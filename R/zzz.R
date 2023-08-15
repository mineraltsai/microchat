"mergeigraph" <- function(g) {
  require(igraph)
  edge_table<-list()
  edge_names<-list()
  node_table<-list()

  for (i in 1:(length(g)-1)) {
    E(g[[i]])$correlation <- E(g[[i]])$weight
    E(g[[i]])$weight <- abs(E(g[[i]])$weight)
    node_table[[i]] = as.data.frame(vertex.attributes(g[[i]]))
    node_table[[i]]$degree = igraph::degree(g[[i]])

    edge_table[[i]] = as.data.frame(edge.attributes(g[[i]]))
    edge_names[[i]] = unlist(strsplit(attr(E(g[[i]]), "vnames"), "|", fixed=T))
    edge_table[[i]]$node1 = edge_names[[i]][seq(from = 1, to = (length(edge_names[[i]]) - 1), by = 2)]
    edge_table[[i]]$node2 = edge_names[[i]][seq(from = 2, to = length(edge_names[[i]]), by = 2)]
    edge_table[[i]]$cor[which(edge_table[[i]]$correlation < 0)] = "-1"
    edge_table[[i]]$cor[which(edge_table[[i]]$correlation > 0)] = "1"

    g[[i+1]]<-g[[i+1]] %>%
      add_vertices(length(node_table[[i]]$name),
                   attr = list(name=node_table[[i]]$name)) %>%
      add_edges(edge_names[[i]],
                attr = list(weight=edge_table[[i]]$correlation))
  }

  gx<- g[[length(g)]]
  E(gx)$correlation <- E(gx)$weight
  E(gx)$weight <- abs(E(gx)$weight)
  edge_table = as.data.frame(edge.attributes(gx))
  edge_names = unlist(strsplit(attr(E(gx), "vnames"), "|", fixed=T))
  edge_table$node1 = edge_names[seq(from = 1, to = (length(edge_names) - 1), by = 2)]
  edge_table$node2 = edge_names[seq(from = 2, to = length(edge_names), by = 2)]

  edge_table$edge<-paste(edge_table$node1,edge_table$node2,sep = "-")
  edge_table<-edge_table[-which(duplicated(edge_table$edge)==TRUE),]

  edge_table<-subset(edge_table, select = c(node1,node2,correlation))
  colnames(edge_table)[3]<-"weight"
  g1<-graph_from_data_frame(edge_table,directed = FALSE, vertices = NULL)

  return(g1)
}

"get.moduleEigengenes" <- function(g1,submchat) {
  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)
  wtc <- igraph::cluster_fast_greedy(g1,NA)
  V(g1)$module<-igraph::membership(wtc)
  node_table = as.data.frame(vertex.attributes(g1))
  aaaa<-t(submchat$otu_table)%>%data.frame()
  MEs0 = WGCNA::moduleEigengenes(aaaa[,V(g1)$name], colors = wtc$membership)$eigengenes
  mes1 <- MEs0
  return(mes1)
}

"microchatEGHeatmap"<-function(mes1,
                               axis.x.rev=FALSE,
                               color_sig=sigcolor,
                               method = "spearman") {
  corr_me_env<-psych::corr.test(mes1,mes1,adjust="fdr",method = method)
  resr<-corr_me_env$r%>%data.frame()
  resp<-corr_me_env$p%>%data.frame()

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

  if (!axis.x.rev) com$variable<-factor(com$variable,levels=(unique(com$variable)))
  if (axis.x.rev) com$variable<-factor(com$variable,levels=rev(unique(com$variable)))
  com$group<-factor(com$group,levels=unique(com$variable))
  com$p.value<-round(com$p.value,2)
  com$sig<-ifelse(com$p.value<0.001,"***",
                  ifelse(com$p.value<=0.01,"**",
                         ifelse(com$p.value<0.05,"*","")))
  com$r[which(com$variable==com$group)]=0
  com$sig[which(com$variable==com$group)]=""

  p1<-ggplot(com, mapping = aes(x=group,y=variable)) +
    geom_tile(linetype="dashed",
              fill="white",aes(color=r))+
    geom_tile(aes(height=r,width=r,fill = r))+
    scale_fill_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    scale_color_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    geom_text(aes(label=sig),vjust=0.75,color="black",angle=45)+
    theme(text = element_text(family = "serif"),
          title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          #panel.background = element_rect(fill = "white",color="white"),
          panel.grid = element_line(color=NA),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(angle=270,vjust = 0.5),
          aspect.ratio = length(unique(com$group))/length(unique(com$variable)))+
    coord_flip()

  return(p1)
}

"microchatEGMatHeatmap"<-function(mes1,p.mat=NULL,
                               axis.x.rev=FALSE,
                               color_sig=sigcolor,
                               method = "spearman") {
  resr<-mes1%>%data.frame()
  resr1<-resr
  resr1$group<-rownames(resr1)
  cxcr<-reshape2::melt(resr1)
  colnames(cxcr)[3]<-"r"

  if (!is.null(p.mat)) {
    resp<-p.mat%>%data.frame()
    resp1<-resp
    resp1$group<-rownames(resp1)
    cxcp<-reshape2::melt(resp1)
    colnames(cxcp)[3]<-"p.value"
    com<-cbind(cxcr,cxcp)
    com<-subset(com,select = c(2,1,3,6))
  } else {
    com<-cxcr
    com<-subset(com,select = c(2,1,3))
    }

  if (!axis.x.rev) com$variable<-factor(com$variable,levels=(unique(com$variable)))
  if (axis.x.rev) com$variable<-factor(com$variable,levels=rev(unique(com$variable)))
  #com$group<-factor(com$group,levels=unique(com$variable))
  if (!is.null(p.mat)) com$p.value<-round(com$p.value,2)
  if (!is.null(p.mat)) com$sig<-ifelse(com$p.value<0.001,"***",
                  ifelse(com$p.value<=0.01,"**",
                         ifelse(com$p.value<0.05,"*","")))
  if (!is.null(p.mat)) com$r[which(com$variable==com$group)]=0
  if (!is.null(p.mat)) com$sig[which(com$variable==com$group)]=""

  if (!is.null(p.mat)) p1<-ggplot(com, mapping = aes(x=group,y=variable)) +
    geom_tile(linetype="dashed",
              fill="white",aes(color=r))+
    geom_tile(aes(height=r,width=r,fill = r))+
    scale_fill_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    scale_color_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    geom_text(aes(label=sig),vjust=0.75,color="grey75",angle=45)+
    theme(text = element_text(family = "serif"),
          title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          #panel.background = element_rect(fill = "white",color="white"),
          panel.grid = element_line(color=NA),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(angle=270,vjust = 0.5),
          aspect.ratio = length(unique(com$group))/length(unique(com$variable)))+
    coord_flip()

  if (is.null(p.mat)) p1<-ggplot(com, mapping = aes(x=group,y=variable)) +
    geom_tile(linetype="dashed",
              fill="white",aes(color=r))+
    geom_tile(aes(height=r,width=r,fill = r))+
    scale_fill_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    scale_color_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    theme(text = element_text(family = "serif"),
          title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          #panel.background = element_rect(fill = "white",color="white"),
          panel.grid = element_line(color=NA),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(angle=270,vjust = 0.5),
          aspect.ratio = length(unique(com$group))/length(unique(com$variable)))+
    coord_flip()

  return(p1)
}


"microchatEGParmHeatmap"<-function(rmnet_param_mult,mes1,
                                   add.anno=TRUE,
                                   add.anno.range=0,
                                   add.anno.alpha=1,
                                   add.sec.param=TRUE,
                                   mytheme=NULL,
                                   sig.text.size=2,
                                   axis.x.rev=FALSE,
                                   color_sig=sigcolor,
                                   method = "spearman",
                                   export_path) {

  export_path<-paste(export_path,"/microbial parameter mining analysis",sep = "")
  dir.create(export_path, recursive = TRUE)

  mat<-rmnet_param_mult[[1]]
  id.sel<-rmnet_param_mult[[2]]
  if (ncol(id.sel)==3) id.sel=id.sel else id.sel[,"trait"]<-"otu"
  id.sel$trait<-factor(id.sel$trait,levels = unique(id.sel$trait))
  id.sel<-id.sel[order(id.sel$trait,id.sel$module),]
  id.sel$me<-paste(id.sel$trait,id.sel$module,sep = "")

  corr_me_env<-psych::corr.test(mat,mes1,adjust="fdr",method = method)
  resr<-corr_me_env$r%>%data.frame()
  resp<-corr_me_env$p%>%data.frame()

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
  #if (length(unique(com$group))==1) com$group<-paste0(id.sel$trait,id.sel$module)[1]

  if (!axis.x.rev) com$variable<-factor(com$variable,levels=(unique(com$variable)))
  if (axis.x.rev) com$variable<-factor(com$variable,levels=rev(unique(com$variable)))
  #com$group<-factor(com$group,levels=unique(com$variable))
  com$p.value<-round(com$p.value,2)
  com$sig<-ifelse(com$p.value<0.001,"***",
                  ifelse(com$p.value<=0.01,"**",
                         ifelse(com$p.value<0.05,"*","")))

  com$remark<-ifelse(com$r>0 & com$sig!="","up",
                     ifelse(com$r<0 & com$sig!="","down","not"))
  p1<-ggplot(com, mapping = aes(x=variable,y=group)) +
    geom_tile(linetype="dashed",
              fill="white",aes(color=r))+
    scale_color_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    geom_tile(aes(height=r,width=r,fill = r))+
    scale_fill_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    geom_text(aes(label=sig),vjust=0.75,color="black",angle=45,size=sig.text.size)+
    theme(text = element_text(family = "serif"),
          title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          #panel.background = element_rect(fill = "white",color="white"),
          panel.grid = element_line(color=NA),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(angle=270,vjust = 0.5),
          aspect.ratio = length(unique(com$group))/length(unique(com$variable)))

  if (add.anno) {
    rowanno<-abs(resr)%>%rowSums()
    colanno<-abs(resr)%>%colSums()

    x.baselen<-length(unique(com$variable))
    rowanno<-rowanno/max(rowanno)


    if (!is.null(names(rowanno))) {
      data.f<-strsplit(names(rowanno),split="(?<=[a-zA-Z])\\s*(?=[0-9])",perl=TRUE)%>%data.frame()
      data.fx<-data.f%>%t()%>%data.frame()
    }

    color.usex<-rep(color_group,2)
    color.usex<-color.usex[1:length(unique(data.fx$X1))]
    names(color.usex)<-unique(data.fx$X1)
    color.use.sel<-color.usex[match(data.fx$X1,names(color.usex))]

    for (i in 1:length(rowanno)) {
      p1<-p1+annotate(geom = "rect",
                      xmin = x.baselen+1-0.4,
                      xmax = x.baselen+1+add.anno.range+rowanno[i]-0.4 ,
                      ymin = i-0.4,
                      ymax = i+0.4,
                      fill=color.use.sel[i],
                      alpha = add.anno.alpha)
    }

    y.baselen<-length(unique(com$group))
    colanno<-colanno/max(colanno)
    for (j in 1:length(colanno)) {
      p1<-p1+annotate(geom = "rect",
                      xmin = j-0.4,
                      xmax = j+0.4 ,
                      ymin = y.baselen+1-0.4,
                      ymax = y.baselen+1+add.anno.range+colanno[j]-0.4,
                      fill=color_group[j],
                      alpha = add.anno.alpha)
    }


    id.sel$pos<-seq(0.5,y.baselen+0.5,length.out=nrow(id.sel))
    neg.or<-1:length(unique(id.sel$me))
    names(neg.or)<- unique(id.sel$me)
    id.sel$neg<-neg.or[match(id.sel$me,names(neg.or))]
    id.sel$type<-grepl('gene',id.sel$trait)
    id.sel$type1<-id.sel$type
    if (length(which(id.sel$type=="TRUE"))!=0) id.sel[which(id.sel$type=="TRUE"),]$type1="italic"
    if (length(which(id.sel$type!="TRUE"))!=0) id.sel[which(id.sel$type!="TRUE"),]$type1="plain"

    if (add.sec.param) for (k in 1:nrow(id.sel)) {

      p1<-p1+annotate(geom = "text",
                 family="serif",
                 x = -3,
                 y =id.sel$pos[k],
                 label = "") +
        annotate(geom = "segment",
                 x = -0.85,
                 xend = -2.8,
                 linetype="dashed",
                 size=0.5,
                 y = id.sel$pos[k],
                 yend = id.sel$neg[k],
                 colour = "grey50")+
        annotate(geom = "text",
                 family="serif",
                 size=2,
                 fontface=id.sel$type1[k],
                 x = -0.4,
                 y =id.sel$pos[k],
                 label = id.sel$name[k])

    }

p1<-p1+annotate(geom = "text",
                family="serif",
            x = x.baselen+1+add.anno.range+max(rowanno)+1,
            y = y.baselen/2,
            vjust=2,
            angle=270,
            label = "Accumulative absolute parametric weight")

p1<-p1+annotate(geom = "text",
                family="serif",
                vjust=2,
            x = x.baselen/2,
            y = y.baselen+1+add.anno.range+max(colanno)+1,
            label = "Accumulative absolute weight")

if (add.sec.param)p1<-p1+
      theme(aspect.ratio = (y.baselen+1+add.anno.range+max(colanno)+1)/(x.baselen+1+add.anno.range+1+max(rowanno)+3))
#p1$theme$aspect.ratio
    if (!add.sec.param) p1<-p1+
      theme(aspect.ratio = (y.baselen+1+add.anno.range+1+max(colanno))/(x.baselen+1+add.anno.range+1+max(rowanno)))
  }
  if (!is.null(mytheme)) p1<-p1+mytheme

  hw.rat<-p1$theme$aspect.ratio

  h.sel<-(y.baselen+1+add.anno.range+1)
  if (h.sel>18) h.sel<-12

  ggsave(paste(export_path,"/microbiota-parameter heatmap (",ncol(mes1)," module eigengene *",length(unique(id.sel$me))," params).pdf",sep = ""),
        height = h.sel,
        width =h.sel/hw.rat,
         p1)
  ggsave(paste(export_path,"/microbiota-parameter heatmap (",ncol(mes1)," module eigengene *",length(unique(id.sel$me))," params).tiff",sep = ""),
         height = h.sel,
         width =h.sel/hw.rat,
         p1)
  cat("图片已保存至相对路径","/",export_path,"","下",sep = "")
  print(p1)
  return(p1)
}


"calcTaxFunCor" <- function(submchat,mes1,
                            min.env.prop =30,
                            ko1.select="Metabolism",
                            method = "spearman",
                            microchatTaxFun) {

  env<-microchatTaxFun$funpred
  envx1<-subset(env,select = c(2:(ncol(submchat$otu_table))))*100

  env<-env[-which(rowSums(envx1)<min.env.prop),]
  env<-subset(env,select = c(2:(ncol(submchat$otu_table)+1),KO1,KO2,KO3))

  if (ko1.select!="all") env<-env[which(env$KO1==ko1.select),] else env<-env

  env<-subset(env,select = c(1:(ncol(submchat$otu_table)),KO3))
  env$KO3<-factor(env$KO3,levels = unique(env$KO3))

  all.index<-group_generate(submchat$otu_table)$sample
  summ_fcbv.all<-data.frame()
  index<-all.index[2]
  for (index in all.index) {
    # calculate the mean and standard error of each group according to the alpha index
    if (index==all.index[1]) {
      summ_fcbv<- env %>%
        group_by(KO3) %>%
        plyr::summarise(
          index=aggregate(as.formula(paste(index, " ~ KO3",sep = "")),
                          data=env,FUN = "sum")) %>%data.frame()
      summ_fcbv<-summ_fcbv$index
      colnames(summ_fcbv)<-c("ko3",index)
      summ_fcbv.all<-rbind(summ_fcbv.all,summ_fcbv)
    } else {
      summ_fcbv<- env %>%
        group_by(KO3) %>%
        plyr::summarise(
          index=aggregate(as.formula(paste(index, " ~ KO3",sep = "")),
                          data=env,FUN = "sum")) %>%data.frame()
      summ_fcbv<-summ_fcbv$index
      colnames(summ_fcbv)<-c("ko3",index)
      summ_fcbv<-summ_fcbv[,-1]
      summ_fcbv.all<-cbind(summ_fcbv.all,summ_fcbv)
    }
  }
  colnames(summ_fcbv.all)[2:ncol(summ_fcbv.all)]<-all.index
  envx<-summ_fcbv.all%>%tibble::column_to_rownames(var = "ko3")%>%
    t()%>%data.frame()


  MEs0<-mes1
  corr_me_env<-psych::corr.test(envx,MEs0,adjust="fdr",method = method)
  resr<-corr_me_env$r%>%data.frame()
  resp<-corr_me_env$p%>%data.frame()

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
  com$group<-rep(summ_fcbv.all$ko3,nrow(com)/length(summ_fcbv.all$ko3))
  simr<-as_adjacency_matrix(g1, type = "both",attr="weight")%>%as.matrix()
  return(list(com=com,corr_me_env=corr_me_env,simr=simr,envx=envx))
}

"calcmicrochatCor" <- function(submchat,mes1,microchatTaxFun,g1,
                               min.env.prop =30,
                               ko1.select="Metabolism",
                               method = "spearman") {

  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)
  wtc <- igraph::cluster_fast_greedy(g1,NA)
  V(g1)$module<-igraph::membership(wtc)
  node_table = as.data.frame(vertex.attributes(g1))

  com<-calcTaxFunCor(submchat,mes1,
                     min.env.prop =min.env.prop,
                     ko1.select=ko1.select,
                     method = method,
                     microchatTaxFun)$com

  corr_me_env<-calcTaxFunCor(submchat,mes1,
                             min.env.prop =min.env.prop,
                             ko1.select=ko1.select,
                             method = method,
                             microchatTaxFun)$corr_me_env

  simr<-calcTaxFunCor(submchat,mes1,
                      min.env.prop =min.env.prop,
                      ko1.select=ko1.select,
                      method = method,
                      microchatTaxFun)$simr


  envx<-calcTaxFunCor(submchat,mes1,
                      min.env.prop =min.env.prop,
                      ko1.select=ko1.select,
                      method = method,
                      microchatTaxFun)$envx

  com$variable<-factor(com$variable,levels=(unique(com$variable)))
  com$p.value<-round(com$p.value,2)
  com$sig<-ifelse(com$p.value<0.001,"***",
                  ifelse(com$p.value<=0.01,"**",
                         ifelse(com$p.value<0.05,"*","")))
  return(list(com=com,
              corr_me_env=corr_me_env,
              simr=simr,
              mes1=mes1,
              envx=envx,
              node_table=node_table,
              submchat=submchat))

}

"microchatEGTaxFunHeatmap" <- function(mcCor,sig.size=5,
                                      color_sig=sigcolor) {
  print(unique(mcCor$com$group))
  p1<-ggplot(mcCor$com, mapping = aes(x=group,y=variable)) +
    geom_tile(linetype="dashed",
              fill="white",aes(color=r))+
    geom_tile(aes(height=r,width=r,fill = r))+
    scale_fill_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    scale_color_gradient2(
      limits=c(-1,1),
      midpoint = 0,
      low = color_sig[1],
      mid = "white",
      high = color_sig[3],
      space = "Lab" ,n.breaks=6)+
    geom_text(aes(label=sig),vjust=0.75,color="grey75",angle=45,size=sig.size)+
    theme(text = element_text(family = "serif"),
          title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          #panel.background = element_rect(fill = "white",color="white"),
          panel.grid = element_line(color=NA),
          #axis.text.y = element_blank(),
          axis.text.x = element_text(angle=270,vjust = 0.5),
          aspect.ratio = length(unique(mcCor$com$group))/length(unique(mcCor$com$variable)))+
    coord_flip()
  return(p1)
}

"calcmicrochatBiomarker" <- function(mcCor,
                                     cor.thres=0.2,
                                     p.thres=0.05,
                                     mod.sel="ME1",
                                     param.sel="Nitrogen metabolism") {

  sel.module<-gsub("ME","",mod.sel)%>%as.numeric()
  gene_module <- mcCor$node_table
  gene_module_select <- subset(gene_module, module == sel.module)$name
  message("\n","OTU belonging to selected module (",mod.sel,") :")
  print(paste(mod.sel,": ",gene_module_select,sep = ""),sep = "")

  #otu matrix according to selected module
  gene_select <- t(t(mcCor$submchat$otu_table)[,gene_module_select])

  #otu co-expression similarity matrix according to selected module
  tom_select <- mcCor$simr[gene_module_select,gene_module_select]

  #otu name according to selected module
  gene_black <- t(mcCor$submchat$otu_table)[ ,gene_module_select]
  identical(gene_black,t(gene_select))

  ###calculate module membership of all otus in the selected module
  ###calculate correlation between otu abundance and selected module eigengene
  geneModuleMembership <- signedKME(gene_black,
                                    mcCor$mes1[mod.sel],
                                    outputColumnName = 'MM')
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),
                                             nrow(mcCor$mes1)))

  ###calculate correlation between otu abundance and selected parameters
  act.name<-mcCor$com$group%>%unique()%>%
    tolower()
  param.sel<-param.sel%>%tolower()
  matnum<-match(param.sel,act.name)
  cat("Selected parameters was ")
  message((mcCor$com$group%>%unique())[matnum])
  geneTraitSignificance <- as.data.frame(cor(gene_black,
                                             mcCor$envx[matnum],
                                             use = 'p'))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),
                                             nrow(mcCor$envx)))

  geneModuleMembership[abs(geneModuleMembership)<cor.thres | MMPvalue>p.thres] <- 0
  geneTraitSignificance[abs(geneTraitSignificance)<cor.thres | MMPvalue>p.thres] <- 0

  select <- cbind(geneModuleMembership, geneTraitSignificance)
  select <- subset(select, geneModuleMembership>=cor.thres & geneTraitSignificance>=cor.thres)
  head(select)
  if (nrow(select)==0) warning("No candidate biomarkers could be found !!!")

  message("Totally finding ",nrow(select), " candidate biomarkers originated from ",mod.sel, " and ",unique(mcCor$com$group)[matnum])
  return(list(otu_select=gene_select,
              mod.sel=mod.sel,
              mat=tom_select,
              select=select,
              selected.otu=rownames(gene_select)))
}


"plotmicrochatBiomarkerBoxplot"<-function(submchat,
                                          sel.otu,
                                          xlabname=xlabname,
                                          ylim.fold=1, ##值越大，y轴显示范围越大
                                          strictmod=TRUE,
                                          method="t.test",
                                          my_comparisons=my_comparisons,
                                          color_group=color_group,
                                          export_path=export_path) {

  export_path<-paste(export_path,"/microbial parameter mining analysis/whole network biomarkers sink/",sep = "")
  dir.create(export_path, recursive = TRUE)


  otu.sel<-submchat$otu_table%>%t()%>%data.frame()
  otu.sel<-otu.sel[,sel.otu]
  params<-otu.sel
  mchat1<-setParamchat(NULL,NULL,NULL,params=params)
  tidymchat1<-tidyMicrochat(mchat1,
                            group.keep=NULL,  ### c("CT","SS","SA","SF","SC)
                            group.rm=NULL,    ### ct
                            sample.keep=NULL, ### c("ct1","ml1")
                            sample.rm=NULL,   ### c("ct1","ml1")
                            taxon.keep=NULL,  ### c("p__firm","c__bacter")
                            taxon.rm=NULL,    ### c("firm","cyan")
                            filter_pollution=TRUE,
                            filter_Archaea=TRUE)
  submchat1<-subsampleMicrochat(tidymchat1,sample.size=NULL,export_path=export_path)
  microchatParamobj<-calcMicrochatParam(submchat1,
                                        export_path=export_path)

  p<-plotMicrochatParamMutiBoxplot(microchatParamobj,
                                   xlabname=xlabname,
                                   geom.line=FALSE,  ##是否添加线条图层
                                   geom.line.size=0.5, ##线条粗细
                                   ylim.fold=ylim.fold, ##值越大，y轴显示范围越大
                                   yaxis.italic=FALSE,
                                   strictmod=strictmod,
                                   method=method,
                                   comparison=my_comparisons,
                                   color_group=color_group,
                                   export_path=export_path)

  return(p)
}

"plotmicrochatBiomarkerBarplot"<-function(submchat,
                                          sel.otu,
                                          xlabname=xlabname,
                                          ylim.fold=1, ##值越大，y轴显示范围越大
                                          strictmod=TRUE,
                                          method="t.test",
                                          my_comparisons=my_comparisons,
                                          color_group=color_group,
                                          export_path=export_path) {

  export_path<-paste(export_path,"/microbial parameter mining analysis/whole network biomarkers sink/",sep = "")
  dir.create(export_path, recursive = TRUE)


  otu.sel<-submchat$otu_table%>%t()%>%data.frame()
  otu.sel<-otu.sel[,sel.otu]
  params<-otu.sel
  mchat1<-setParamchat(NULL,NULL,NULL,params=params)
  tidymchat1<-tidyMicrochat(mchat1,
                            group.keep=NULL,  ### c("CT","SS","SA","SF","SC)
                            group.rm=NULL,    ### ct
                            sample.keep=NULL, ### c("ct1","ml1")
                            sample.rm=NULL,   ### c("ct1","ml1")
                            taxon.keep=NULL,  ### c("p__firm","c__bacter")
                            taxon.rm=NULL,    ### c("firm","cyan")
                            filter_pollution=TRUE,
                            filter_Archaea=TRUE)
  submchat1<-subsampleMicrochat(tidymchat1,sample.size=NULL,export_path=export_path)
  microchatParamobj<-calcMicrochatParam(submchat1,
                                        export_path=export_path)

  p<-plotMicrochatParamMutiBarplot(microchatParamobj,
                                   xlabname=xlabname,
                                   ylim.fold=ylim.fold, ##值越大，y轴显示范围越大
                                   yaxis.italic=FALSE,
                                   strictmod=strictmod,
                                   method=method,
                                   comparison=my_comparisons,
                                   color_group=color_group,
                                   export_path=export_path)


  return(p)
}

"rmnet_param" <- function(param,cor.method = "spearman") {
  otu_rare<-param
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

  if (as.numeric(matchnum) != as.numeric(allnum)) split_otu <- lapply(
    sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),
    FUN = function(x){otu_rare[,x]})

  if (as.numeric(matchnum) == as.numeric(allnum)) split_otu <- apply(
    sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),2,
    FUN = function(x){otu_rare[,x]})

  otu<-otu_filter(split_otu,filter_num = 1,filter = FALSE,self = FALSE)


  mat <- function(g) {
    gmat<-lapply(g, function(x){
      gg<-as.matrix(x)
      return(gg)
    })
    return(gmat)
  }
  ###subgroup
  otus<-mat(otu)

  nor <- lapply(otus,function(x){
    new<-scale(t(x))##
    occor.r = cor(new,method = cor.method,use = "pairwise.complete.obs")
    diag(occor.r) <- 0
    occor.r[is.na(occor.r)]=0
    return(occor.r)
  })

  net <- lapply(nor,function(x){
    cleaned.matrix <- x
    g = graph_from_adjacency_matrix(cleaned.matrix ,mode="undirected",weighted=TRUE,diag=FALSE)
    result<-0
    name<-names(x)
    return(list(g=g,cleaned.matrix=cleaned.matrix,result=0,name=name))
  })
  g<-lapply(net, function(x){
    gg<-x$g
    return(gg)
  })
  return(list(net=net,g=g))
}

"rmnet_micro" <- function(g,abun,cor.method = "spearman") {
  colname<-colnames(abun)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)

  group_num<-length(trt_id)

  groupname<-substr(colnames(abun),start = 1,stop = 2)%>%unique()
  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(abun))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()

  if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otu <- lapply(
      apply(
        sapply(trt_id,function(x){grep(x,colnames(abun))}),2,
        FUN = function(x){abun[,x]}),
      function(x){x[-(which(rowSums(x)==0)),]})
  } else {
    split_otu <- lapply(
      lapply(
        sapply(trt_id,function(x){
          grep(x,colnames(abun))
        }),FUN = function(x){
          abun[,x]
        }),
      function(x){x[-(which(rowSums(x)==0)),]})
  }

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

  return(eige)
}

"rmnet_param_eigene" <- function(param,cor.method) {

  #param<-submchat$param_table[[1]]
  param<-param%>%t()%>%data.frame()
  rmnet_xparam<-rmnet_param(param,cor.method = cor.method)
  abun<-param
  g<-rmnet_xparam$g

  colname<-colnames(abun)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)

  split_otu <- apply(
    sapply(trt_id,function(x){grep(x,colnames(abun))}),2,
    FUN = function(x){abun[,x]})

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

  return(eige)

}

"get.paramat.eigen"<-function(submchat,method) {
  paramx<-submchat$param_table
  if (class(paramx)=="list") {
    rmnet_param_mutieigene<-data.frame()
    #tk=10
    for (tk in 1:length(paramx)) {
      if (tk==1) {
        trats<-names(paramx)[[1]]
        param<-paramx[[1]]
        rmnet_param_xeigene<-rmnet_param_eigene(param,method)
        colnames(rmnet_param_xeigene)<-paste(trats,1:ncol(rmnet_param_xeigene),sep = "")
        rmnet_param_mutieigene<-rbind(rmnet_param_mutieigene,rmnet_param_xeigene)
      } else {
        trats<-names(paramx)[[tk]]
        param<-paramx[[tk]]
        rmnet_param_xeigene<-rmnet_param_eigene(param,method)
        colnames(rmnet_param_xeigene)<-paste(trats,1:ncol(rmnet_param_xeigene),sep = "")

        {
          rmnet_param_mutieigene<-cbind(rmnet_param_mutieigene,rmnet_param_xeigene)
        }
      }

    }
  } else {
    param<-paramx
    rmnet_param_xeigene<-rmnet_param_eigene(param,method)
    colnames(rmnet_param_xeigene)<-paste("param",1:ncol(rmnet_param_xeigene),sep = "")
    rmnet_param_mutieigene<-rmnet_param_xeigene
  }
  message("selected param table was ",class(paramx))
  return(rmnet_param_mutieigene)
}

"rmnet_param.eigene" <- function(param,cor.method) {

  #param<-submchat$param_table[[1]]
  param<-param%>%t()%>%data.frame()
  rmnet_xparam<-rmnet_param(param,cor.method = cor.method)

  gx<-mergeigraph(rmnet_xparam$g)
  E(gx)$correlation <- E(gx)$weight
  E(gx)$weight <- abs(E(gx)$weight)
  wtc <- igraph::cluster_fast_greedy(gx,NA)
  V(gx)$module<-igraph::membership(wtc)
  node_table = as.data.frame(vertex.attributes(gx))
  mesx<-get.paramEigengenes(gx,param)

  return(list(mes=mesx,node_table=node_table,g=gx))

}

"get.multparam.eigen" <- function(submchat,method) {
  paramx<-submchat$param_table
  if (class(paramx)=="list") {
    rmnet_param_mutieigene<-data.frame()
    rmnet_param_multnt<-data.frame()
    #tk=1
    for (tk in 1:length(paramx)) {
      if (tk==1) {
        trats<-names(paramx)[[1]]
        param<-paramx[[1]]
        rmnet_param_xeigene<-rmnet_param.eigene(param,cor.method=method)$mes
        rmnet_param_nt<-rmnet_param.eigene(param,cor.method=method)$node_table
        rmnet_param_nt$trait<-trats
        rmnet_param_multnt<-rbind(rmnet_param_multnt,rmnet_param_nt)

        colnames(rmnet_param_xeigene)<-paste(trats,1:ncol(rmnet_param_xeigene),sep = "")
        rmnet_param_mutieigene<-rbind(rmnet_param_mutieigene,rmnet_param_xeigene)

      } else {
        trats<-names(paramx)[[tk]]
        param<-paramx[[tk]]
        rmnet_param_xeigene<-rmnet_param.eigene(param,cor.method=method)$mes
        colnames(rmnet_param_xeigene)<-paste(trats,1:ncol(rmnet_param_xeigene),sep = "")
        rmnet_param_nt<-rmnet_param.eigene(param,cor.method=method)$node_table
        rmnet_param_nt$trait<-trats

        {
          rmnet_param_multnt<-rbind(rmnet_param_multnt,rmnet_param_nt)
          rmnet_param_mutieigene<-cbind(rmnet_param_mutieigene,rmnet_param_xeigene)
        }
      }

    }
  } else {
    param<-paramx
    rmnet_param_xeigene<-rmnet_param.eigene(param,cor.method=method)$mes
    colnames(rmnet_param_xeigene)<-paste("param",1:ncol(rmnet_param_xeigene),sep = "")
    rmnet_param_mutieigene<-rmnet_param_xeigene

    rmnet_param_nt<-rmnet_param.eigene(param,cor.method=method)$node_table
    rmnet_param_nt$trait<-"param"
    rmnet_param_multnt<-rmnet_param_nt

  }
  message("selected param table was ",class(paramx))
  return(list(rmnet_param_mutieigene=rmnet_param_mutieigene,
              rmnet_param_multnt=rmnet_param_multnt))
}

"get.paramEigengenes" <- function(gx,abun) {

  g1<-gx
  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)
  wtc <- igraph::cluster_fast_greedy(g1,NA)
  V(g1)$module<-igraph::membership(wtc)
  node_table = as.data.frame(vertex.attributes(g1))
  aaaa<-abun%>%t()%>%data.frame()
  MEs0 = WGCNA::moduleEigengenes(aaaa[,V(g1)$name], colors = wtc$membership)$eigengenes
  mes1 <- MEs0
  return(mes1)
}

"get.paraModule" <- function(gx) {

  g1<-gx
  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)
  wtc <- igraph::cluster_fast_greedy(g1,NA)
  V(g1)$module<-igraph::membership(wtc)
  node_table = as.data.frame(vertex.attributes(g1))

  return(node_table)
}

"calcEGParmBiomarker" <- function(g1,submchat,
                                  rmnet_param_mult,
                                  cor.thres=0.2,
                                  p.thres=0.05,
                                  otu.mod.sel="ME2",
                                  param.mod.sel="gut_gene1") {
require(WGCNA)
  mod.sel=otu.mod.sel
  param.sel=param.mod.sel

  sel.module<-gsub("ME","",mod.sel)%>%as.numeric()
  gene_module<- get.paraModule(g1)
  gene_module_select <- subset(gene_module, module == sel.module)$name
  message("\n","OTU belonging to selected module (",mod.sel,") :")
  print(paste(mod.sel,": ",gene_module_select,sep = ""),sep = "")

  #otu matrix according to selected module
  gene_select <- t(t(submchat$otu_table)[,gene_module_select])

  #otu co-expression similarity matrix according to selected module
  simr<-as_adjacency_matrix(g1, type = "both",attr="weight")%>%as.matrix()
  tom_select <- simr[gene_module_select,gene_module_select]

  #otu name according to selected module
  gene_black <- t(submchat$otu_table)[ ,gene_module_select]
  identical(gene_black,t(gene_select))

  mes1<-get.moduleEigengenes(g1,submchat)
  ###calculate module membership of all otus in the selected module
  ###calculate correlation between otu abundance and selected module eigengene
  geneModuleMembership <- signedKME(gene_black,
                                    mes1[mod.sel],
                                    outputColumnName = 'MM')
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),
                                             nrow(mes1)))

  ###calculate correlation between otu abundance and selected parameters
  com<-rmnet_param_mult[[1]]
  param.select<-rmnet_param_mult[[2]]
  if (ncol(param.select)==3) param.select<-param.select else param.select[,"trait"]<-"ME"
  param.select$param<-paste(param.select$trait,param.select$module,sep = "")

  param.select<-param.select[which(param.select$param==param.sel),]

  act.name<-colnames(com)%>%unique()%>%
    tolower()
  param.sel<-param.sel%>%tolower()
  matnum<-match(param.sel,act.name)
  cat("Selected parameters was ")
  message(colnames(com)[matnum])
  geneTraitSignificance <- as.data.frame(cor(gene_black,
                                             com[matnum],
                                             use = 'p'))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),
                                             nrow(com)))

  geneModuleMembership[abs(geneModuleMembership)<cor.thres | MMPvalue>p.thres] <- 0
  geneTraitSignificance[abs(geneTraitSignificance)<cor.thres | MMPvalue>p.thres] <- 0

  select <- cbind(geneModuleMembership, geneTraitSignificance)
  select <- subset(select, abs(geneModuleMembership)>=cor.thres & abs(geneTraitSignificance)>=cor.thres)
  head(select)
  if (nrow(select)==0) warning("No candidate biomarkers could be found !!!")

  gene_select.cand<-gene_select[match(rownames(select),rownames(gene_select)),]
  message("MM1 indicated module membership")
  message("Totally finding ",nrow(select), " candidate biomarkers originated from ",mod.sel, " and ",unique(colnames(com))[matnum])
  message("Totally finding ",nrow(param.select), " candidate parametric indice(s) originated from ",param.sel)
  return(list(otu_select.all=gene_select,
              otu_select=gene_select.cand,
              param_select=param.select,
              mod.sel=mod.sel,
              param.sel=param.sel,
              sim.matrix=tom_select,
              candidate.otu.eigen=select,
              selected.param=param.select$name,
              selected.otu.all=rownames(gene_select),
              selected.otu=rownames(select)))
}


"get.param.eigen" <- function(rmnet_param_mult,param.mod.sel="gut_gene1") {
  dat.f <-rmnet_param_mult[[1]][,param.mod.sel]%>%data.frame()
  colnames(dat.f)<-param.mod.sel
  rownames(dat.f)<-rownames(rmnet_param_mult[[1]])
  rmnet_param_mult[[1]]<-dat.f
  rmnet_param_mult[[2]]<-rmnet_param_mult[[2]][which(paste0(rmnet_param_mult[[2]]$trait,rmnet_param_mult[[2]]$module)%in%param.mod.sel),]
  return(rmnet_param_mult)
}



"get.param.eigen.t" <- function(submchat,
                                rmnet_param_multx) {
  dafca<-submchat$param_table[which(names(submchat$param_table)==unique(rmnet_param_multx$rmnet_param_multnt$trait))]
  data.n<-rmnet_param_multx$rmnet_param_multnt

  data.selx<-data.frame()
  for (j in 1:length(unique(data.n$trait))) {
    if (j==1) {
      param.sel<-data.n[which(unique(data.n$trait)[j]==data.n$trait),]$name
      data.sel<-dafca[[which(names(dafca)==unique(data.n$trait)[j])]]
      data.sel<-data.sel[,param.sel]
      data.selx<-rbind(data.selx,data.sel)
    }
    if (j!=1) {
      param.sel<-data.n[which(unique(data.n$trait)[j]==data.n$trait),]$name
      data.sel<-dafca[[which(names(dafca)==unique(data.n$trait)[j])]]
      data.sel<-data.sel[,param.sel]
      {
        data.selx<-cbind(data.selx,data.sel)
      }
    }
  }
  data.selw<-list()
  data.selw[[1]]<-data.selx
  data.selw[[2]]<-rmnet_param_multx$rmnet_param_multnt
  colnames(data.selw[[1]])[duplicated(colnames(data.selw[[1]]))]<-paste0(colnames(data.selw[[1]])[duplicated(colnames(data.selw[[1]]))],
                                                                         length(duplicated(colnames(data.selw[[1]]))))
  return(data.selw)
}


"get.lefse.resName" <- function(microchatLefseobj) {
  microchatLefseobj<-microchatLefseobj%>%as.list()
  asda<-microchatLefseobj$res[which(substr(microchatLefseobj$res$f,start = 1,stop = 3)=="s__"),]
  wc<-asda[which(substr(asda$f,start = 1,stop = 9)=="s__un_g__"),]
  otu3.n<-taxon[which(paste0(taxon$Species,"s__un_",taxon$Genus) %in% wc$f),]%>%rownames()
  asda<-asda[which(substr(asda$f,start = 1,stop = 5)!="s__un"),]
  otu1.n<-taxon[which(taxon$Species %in% asda$f),]%>%rownames()
  asdx<-microchatLefseobj$res[which(substr(microchatLefseobj$res$f,start = 1,stop = 3)=="OTU"),]
  otu2.n<-asdx$f%>%as.character()

  diff.dat<-c(otu1.n,otu2.n,otu3.n)
  return(diff.dat)
}


"calcMicrochatRF" <- function(submchat) {
  testData<-submchat$otu_table
  testData<-testData%>%t()%>%data.frame()
  testData <- testData %>%
    mutate(treatments = substr(rownames(testData),start = 1,stop = 2))

  testData2 <- testData %>%
    mutate(treatments = as.factor(treatments))

  set.seed(123)
  otu_rf= randomForest::randomForest(treatments ~ ., data = testData2, importance=TRUE, proximity=TRUE)

  imp_otu <- as_tibble(round(importance(otu_rf), 2),rownames = "OTUid") %>%
    arrange(desc(MeanDecreaseAccuracy))

  ncol(testData2)
  myotu= testData2[,-ncol(testData2)]
  set.seed(123)

  result<-randomForest::rfcv(myotu, testData2$treatments, cv.fold=5, scale = "log", step = 0.9)
  result1<-result
  result1$n.var

  plot(result$n.var, result$error.cv, log="x", type="o", lwd=2)
  biomarker.num<-which(result$error.cv==min(result$error.cv))%>%names()%>%as.numeric()%>%min()
  return(list(biomarker.num=biomarker.num,
              result=result,
              result1=result1,myotu=myotu,
              imp_otu=imp_otu,
              testData2=testData2))
}

"get.RFThres" <- function(rf.dat,main_theme) {

  result<-rf.dat$result
  error.cv <- data.frame(num = result$n.var,
                         error.1 =  result$error.cv)
  myotu<-rf.dat$myotu
  testData2<-rf.dat$testData2
  for (i in 124:127){
    print(i)
    set.seed(i)
    result= randomForest::rfcv(myotu, testData2$treatments, cv.fold=5, scale = "log", step = 0.9)
    error.cv = cbind(error.cv, result$error.cv)
  }

  n.var <- error.cv$num
  error.cv <- error.cv[,-1]
  colnames(error.cv)<- paste('err',1:5,sep='.')#添加列名
  err.mean <-apply(error.cv,1,mean)#求5列的平均值
  allerr<-data.frame(num=n.var,err.mean=err.mean,error.cv)

  ggplot() +
    geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black',size=0.5) +
    geom_vline(xintercept = rf.dat$biomarker.num, colour='black', lwd=0.36, linetype=2) +
    coord_trans(x = "log2") +
    scale_x_continuous(breaks = c(1, 5, 10, 20, 30, 50, 100, 300, max(rf.dat$result$n.var))) + # , max(allerr$num)
    labs( x='Number of OTUs ', y='Cross-validation error rate') +
    annotate("text", x = rf.dat$biomarker.num, y = max(allerr$err.mean),
             family="serif",
             label=paste("Optimal = ", rf.dat$biomarker.num, sep=""))+
    main_theme->p

  ###2
  optimal = allerr$num[which(allerr$err.mean==min(allerr$err.mean))]%>%min()
  ggplot() +
    geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey',size=0.5) +
    geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black',size=0.5) +
    geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype=2) +
    coord_trans(x = "log2") +
    scale_x_continuous(breaks = c(1, 5, 10, 20, 30, 50, 100, 300, max(rf.dat$result$n.var))) + # , max(allerr$num)
    labs( x='Number of OTUs ', y='Cross-validation error rate') +
    annotate("text", x = optimal, y = max(allerr$err.mean),
             family="serif",
             label=paste("Optimal = ", optimal, sep=""))+
    main_theme ->p1
  pp<-p+p1
  print(pp)
  imp_otu<-rf.dat$imp_otu
  return(list(p=p,p1=p1,imp_otu=imp_otu,
              optimal.r=optimal,optimal.f=rf.dat$biomarker.num))
}

"plotmicrochatRFbiomarker" <- function(opt.num,main_theme) {
  opt.num.sel<-c(opt.num$optimal.r,opt.num$optimal.f)%>%as.numeric()
  imp_otu<-opt.num$imp_otu
  imp_otu[1:max(opt.num.sel),] %>%
    rstatix::select(OTUid,MeanDecreaseAccuracy) %>%
    arrange(MeanDecreaseAccuracy) %>%
    mutate(OTUid = forcats::fct_inorder(OTUid)) %>%
    ggplot(aes(x = OTUid, y = MeanDecreaseAccuracy))+
    geom_bar(stat = "identity",width = 0.75)+
    labs(x = "", y = "Mean decrease accuracy")+
    coord_flip()+
    main_theme+
    theme(axis.text.y = element_text(size = 8))->p2

  imp_otu[1:min(opt.num.sel),] %>%
    rstatix::select(OTUid,MeanDecreaseAccuracy) %>%
    arrange(MeanDecreaseAccuracy) %>%
    mutate(OTUid = forcats::fct_inorder(OTUid)) %>%
    ggplot(aes(x = OTUid, y = MeanDecreaseAccuracy))+
    geom_bar(stat = "identity",width = 0.75)+
    labs(x = "", y = "Mean decrease accuracy")+
    coord_flip()+
    main_theme+
    theme(axis.text.y = element_text(size = 8))->p3
  print(p2|p3)
  return(list(p.f=p3,p.r=p2))

}
