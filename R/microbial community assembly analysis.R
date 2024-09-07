"calcMicrochatAssemblyNiche" <- function(submchat,
                                    export_path="cs2/microbial assembly analysis") {
  export_path<-paste(export_path,"/data_microbiome/microbial assembly analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,"/Ecological niche width/data",sep = ""), recursive = TRUE)

  if (class(submchat)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  sample.size=submchat$otu_table%>%colSums()
  if (length(unique(sample.size))==1) sample.size<-sample.size[1] else sample.size<-max(sample.size)
  otutab<-submchat$otu_table
  taxtab<-submchat$taxon_table%>%tibble::rownames_to_column(var="name")

  split_otu<-calc_indivGroup_mean(otutab)
  split_otux<-calc_indivGroup_abun(otutab)

  split_otux<-lapply(split_otux, function(x){
    x<-x[-which(rowMeans(x)/sample.size<0.00002),]
  })

  set.seed(1023)
  spec_gen <- lapply(split_otux, function (y) {
    spec.tt<- EcolUtils::spec.gen(t(y),
                                  niche.width.method = 'levins',
                                  perm.method = 'quasiswap',
                                  n = 100,
                                  probs = c(0.025, 0.975))%>%tibble::rownames_to_column(var = "name")
    niche.tt<- spaa::niche.width(t(y), method = "levins")%>%t()%>%data.frame()%>%
      tibble::rownames_to_column(var = "name")
    colnames(niche.tt)[2]<-"niche"
    spec.tt<-merge(spec.tt,niche.tt,by="name")
  })

  thres <- lapply(split_otux, function (y) {
    niche.tt.all<-spaa::niche.width(t(y),
                      method = "levins")%>%t()%>%data.frame()%>%
      tibble::rownames_to_column(var = "name")
    colnames(niche.tt.all)[2]<-"niche"
    niche.outlier.high<-fun.outlier(niche.tt.all$niche)$outlier.high
  })%>%data.frame()

  niche.outlier.high<-rep(min(thres),ncol(thres))

  niche.wid.x <- lapply(split_otux, function (y) {
    niche.tt<- spaa::niche.width(t(y), method = "levins")%>%t()%>%data.frame()%>%
      tibble::rownames_to_column(var = "name")
    colnames(niche.tt)[2]<-"niche"
    niche.tt<-niche.tt
  })

  for (tt in 1:length(spec_gen)) {
    split_otu[[tt]]<-split_otu[[tt]]%>%data.frame()%>%tibble::rownames_to_column(var = "name")
    colnames(split_otu[[tt]])[2]<-"mean"
    spec_gen[[tt]]<-merge(spec_gen[[tt]],split_otu[[tt]],by="name")
    spec_gen[[tt]]<-merge(spec_gen[[tt]],taxtab,by="name")
    if (!is.null(export_path)) {
      file2=paste(export_path,"/Ecological niche width/data/",names(spec_gen)[tt],"_niche width & generalist and specialist.txt",sep = "")
      write.table(spec_gen[[tt]],file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
      cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment has been exported to ","/",export_path,"",sep = "","\n")
    } else {
      cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment was not exported to local path","",sep = "","\n")
    }
  }

  sepec.var<-c("GENERALIST","NON SIGNIFICANT","SPECIALIST")
  spe.gen.prop<-sapply(spec_gen, function(x) {
    prop<-(table(x$sign)/nrow(x))%>%t()%>%data.frame()
    if(nrow(prop)==2) {
      diff.row<-setdiff(c(1,2,3),match(unique(prop$Var2),sepec.var))
      diff.role<-sepec.var[diff.row]
      propx<-rbind(prop,prop[1,])
      propx[match(unique(prop$Var2),sepec.var),]<-prop
      propx.diff<-c("A",diff.role,0)
      propx$Var2<-factor(propx$Var2,levels = sepec.var)
      propx[diff.row,]<-propx.diff
      propx$Var2[diff.row]<-diff.role
      prop<-propx
    } else if (nrow(prop)==1) {
      stop("Only one niche role could be detected!!!")
    }
    prop$Freq<-sprintf("%0.2f",round(prop$Freq*100,2))
    prop$Freq<-paste(prop$Freq,"%",sep = "")
  })%>%data.frame()
  rownames(spe.gen.prop)<-unique(spec_gen[[1]]$sign)[order(unique(spec_gen[[1]]$sign))]
  spe.gen.prop<-spe.gen.prop%>%tibble::rownames_to_column(var = "Metabolism type")

  if (!is.null(export_path)) {
    file2=paste(export_path,"/Ecological niche width/data/Generalist and specialist proportion",".txt",sep = "")
    write.table(spe.gen.prop,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
    cat("\n","Generalist and specialist proportion ","has been exported to ","/",export_path,"",sep = "","\n")
  } else {
    cat("\n","Generalist and specialist proportion ","was not exported to local path","",sep = "","\n")
  }

  pp<-ggpubr::ggtexttable(spe.gen.prop,rows = NULL)%>%
    ggpubr::tab_add_title(text = "Generalist and specialist proportion",
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic") %>%
    ggpubr::tab_add_title(text = "Microbial community assembly analysis",
                          family = "serif",
                          size = 10, face = "bold.italic")%>%
    ggpubr::tab_add_footnote(text = "--Generalist: more flexible in metabolism",
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic") %>%
    ggpubr::tab_add_footnote(text = "--Specialist: less flexible in metabolism",
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic") %>%
    ggpubr::tab_add_footnote(text = "*All data processed by R were presented.",
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic")

  print(pp)
  spec_gen.all<-list(spec_gen=spec_gen,spe.gen.prop=spe.gen.prop,niche.width=niche.wid.x,high.thres=round(niche.outlier.high,0))
  class(spec_gen.all)<-"microchat"
  cat("\n","Based on the niche breadth (B) distribution of OTUs from all treatment, the prescribed minimum in upper outlier area of niche breadth (B) was ",niche.outlier.high[1],",",sep = "")
  message(" which was rounded up to ", round(niche.outlier.high,0)[1],".")
  message("\n","Remove OTUs whose mean relative abundance < 2*10^-5")
  return(spec_gen.all)
}

"plotMicrochatAssemblyNiche"<-function(spec_gen,nrow=1,ncol=3,
                                       xlabname=NULL,
                                       specified.thres=NULL,
                                       sel.shape=c(21,19,23),
                                       color_type=c("red","white","blue"),
                                       nice_layout=FALSE,
                                       color_background=NA,
                                       export_path="cs2/microbial assembly analysis") {
  export_path<-paste(export_path,"/data_microbiome/microbial assembly analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,"/Ecological niche width",sep = ""), recursive = TRUE)

  spec_gen.x<-spec_gen$spec_gen
  spec_gen.prop.x<-spec_gen$spe.gen.prop
  message("specified.thres could be a vector, whose amount is equal with treatment.")
  if (is.null(specified.thres)) high.thres<-spec_gen$high.thres else high.thres<-specified.thres
  pp<-list()
  for (tt in 1:length(spec_gen.x)) {
    sel.g<-names(spec_gen.x)[[tt]]
    if (!is.null(xlabname)) sel.gn<-xlabname[tt] else sel.gn<-sel.g
    sel.prop<-subset(spec_gen.prop.x,select=c(1,match(sel.g,colnames(spec_gen.prop.x))))
    spec_gen.t<-spec_gen.x[[tt]]

    ### high threshold identified by kw.test with the form of mean+std
    spec_gen.t$ntype<-0
    spec_gen.t$ntype<-ifelse(spec_gen.t$niche>=high.thres[tt],"high",
                             ifelse(spec_gen.t$niche<1.5,"low","none"))


    names(color_type)<-c("SPECIALIST","NON SIGNIFICANT","GENERALIST")
    niche.od<-c("SPECIALIST","NON SIGNIFICANT","GENERALIST")
    aaa<-table(spec_gen.t$ntype)
    aaa<-aaa[match(names(aaa),niche.od)]
    bbb<-paste(sprintf("%0.2f",round(aaa/sum(aaa)*100,2)),"%",sep = "")

    names(sel.shape)<- c("SPECIALIST","NON SIGNIFICANT","GENERALIST")

    spe.num<-which(spec_gen.t$sign=="SPECIALIST")%>%length()
    gen.num<-which(spec_gen.t$sign=="GENERALIST")%>%length()
    p<-ggplot(spec_gen.t,aes(mean,niche,color=sign))+
      geom_point(aes(shape=sign),size=1.5)+
      scale_shape_manual(values = sel.shape)+
      scale_color_manual(values = color_type)+
      geom_hline(yintercept=1.5,color=color_type[1],linetype="dashed")+
      geom_hline(yintercept=high.thres[tt],color=color_type[3],linetype="dashed")+
      annotate(geom="text",label=paste("Generalist: ",gen.num," (",sel.prop[1,2],")",sep = ""),
               x=(max(spec_gen.t$mean)-min(spec_gen.t$mean))/2,y=high.thres[tt]+0.3,
               size=3,fontface="bold",
               family="serif",color=color_type[1])+
      annotate(geom="text",label=paste("Specialist: ",spe.num," (",sel.prop[3,2],")",sep = ""),
               x=(max(spec_gen.t$mean)-min(spec_gen.t$mean))/2,y=1.2,
               size=3,fontface="bold",
               family="serif",color=color_type[3])+
      labs(x="Mean abundance",y="Niche breadth (B)",title = sel.gn)+
      theme(panel.background = element_rect(color = color_background,fill = NA),
            text = element_text(family = "serif"),
            aspect.ratio = 1,
            title =  element_text(size = 10,face="bold"),
            axis.title = element_text(size = 12,face="bold"),
            axis.text = element_text(size = 10,face="bold"),
            panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5),
            legend.position = "none")
    if (nice_layout) p<-p+theme(
      panel.border = element_blank()
    )+annotate(geom = "segment",size=2,lineend = "round",
               x = 0, xend =max(spec_gen.t$mean),
               y = -Inf, yend = -Inf)+
      annotate(geom = "segment",size=2,lineend = "round",
               y = min(spec_gen.t$niche), yend = max(spec_gen.t$niche), x = -Inf+1, xend = -Inf+1)+
      annotate(geom = "segment",size=2,lineend = "round",
               x = 0, xend = max(spec_gen.t$mean),
               y = Inf, yend = Inf)+
      annotate(geom = "segment",size=2,lineend = "round",
               y =  min(spec_gen.t$niche), yend = max(spec_gen.t$niche), x = Inf+1, xend = Inf+1)
    if (!is.null(export_path)) {
      file2=paste(export_path,"/Ecological niche width/",sel.g,"_niche breadth.pdf",sep = "")
      ggsave(file2,
             units = "cm",
             width = 21/3,
             height = 21*p$theme$aspect.ratio/3,
             p)
      cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment has been exported to ","/",export_path,"",sep = "","\n")
    } else {
      cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment was not exported to local path","",sep = "","\n")
    }
    pp[[tt]]<-p
  }

  if (!is.null(ncol)){
    width.sel<-ncol
    if ((length(pp)/width.sel)%%1==0) {
      height.sel<-(length(pp)/width.sel)%>%as.integer()
    } else {
      height.sel<-(length(pp)/width.sel)%>%as.integer()+1
    }
  } else {
    width.sel<-ncol_layout(length(pp))[[1]]%>%length()
    height.sel<-ncol_layout(length(pp))%>%length()
  }
  if (!is.null(nrow)) xxxxpx<-nrow
  mp<-patchwork::wrap_plots(pp,ncol = width.sel)

  if (!is.null(export_path)) {
    file2=paste(export_path,"/Ecological niche width/","Muti-group niche breadth.pdf",sep = "")
    rownum<-length(pp)/ncol
    if (rownum>as.integer(rownum)) nrows<-as.integer(rownum)+1 else nrows<-as.integer(rownum)
    ggsave(file2,
           units = "cm",
           width = 7*ncol,height = 7*nrows,
           mp)
    cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment has been exported to ","/",export_path,"",sep = "","\n")
  } else {
    cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment was not exported to local path","",sep = "","\n")
  }

  message("Hollow Circle--Specialist Circle--Non-significant Square--Generalist")
  return(mp)
}

"plotMicrochatNicheWidth" <- function(spec_gen,
                                      color_group=colorCustom(4,pal = "set3"),
                                      panel.border=c("none","normal","all","axis"),
                                      axis.x.angle.adjust=FALSE,
                                      xlabname=NULL,
                                      mytheme=NULL,
                                      color_background=NA,
                                      export_path="cs2/microbial assembly analysis") {
  panel.border<-match.arg(panel.border)
  export_path<-paste(export_path,"/data_microbiome/microbial assembly analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,"/Ecological niche width",sep = ""), recursive = TRUE)

  niche.width<-spec_gen$niche.width
  maxsize<-sapply(niche.width, function(y){max(y$niche[which(y$niche!="Inf")])})%>%data.frame()%>%tibble::rownames_to_column(var = "group")

  xx<-niche.width[[1]]
  colnames(xx)[1]<-"group"
  xx$group<-names(niche.width)[1]
  for (tt in 2:length(niche.width)) {
    colnames(niche.width[[tt]])[1]<-"group"
    niche.width[[tt]]$group<-names(niche.width)[tt]
    xx<-rbind(xx,niche.width[[tt]])
  }
  if (length(which(xx$niche=="Inf"))==0) xx<-xx else xx<-xx[-which(xx$niche=="Inf"),]

  input.data<-reshape2::melt(xx)
  input.data1<-kruskalmuticomp.t(input.data)
  input.data1$group<-factor(input.data1$group,levels = names(niche.width))
  input.data1<-input.data1[order(input.data1$group),]
  input.data1<-merge(input.data1,maxsize,by="group")
  colnames(input.data1)[6]<-"max"

  xx$group<-factor(xx$group,levels = names(niche.width))
  names(color_group)<-names(niche.width)
  p<-ggplot(xx,aes(group,niche))+
    geom_errorbar(data = input.data1,
                  aes(x = group, y = mean, group = group,
                      ymin = mean-std, ymax = mean+std),
                  colour = "black",width=.25)+
    geom_boxplot(outlier.shape = NA,
                 width = 0.5,
                 color = "white")+
    geom_text(data = input.data1,aes(x=group,y=max,label=Letters),
              size=4,vjust=-0.5,family="serif")+
    geom_jitter(data=xx,size=1,alpha=0.75,aes(color=group),
                width = 0.25)+
    geom_point(aes(color=group),size=1)+
    scale_color_manual(values = color_group)+
    scale_fill_manual(values = color_group)+
    labs(y="Niche breadth (B)")+
    theme(
          text = element_text(family = "serif"),aspect.ratio = 1,
          legend.position = "none")

    for (i in 1:nrow(input.data1)) p<-p + annotate(geom="text",
            label=sprintf("%0.2f",round(input.data1$mean[i],2)),
             x=i,y=0,family="serif",color="black",size=3)

  p<-p+
    theme(
      panel.border = element_blank(),
      axis.ticks.length = unit(0.2,"lines"),
      axis.ticks = element_line(color='black'),
      panel.background = element_rect(fill = color_background,color = color_background),
      legend.position = "none",aspect.ratio = 1)

  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-input.data1$group%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p<-p+ scale_x_discrete(labels = xlabname)
  }
  if (panel.border=="all") p<-p+theme(aspect.ratio = 1,
                                      title =  element_text(size = 7,face="bold"),
                                      axis.title = element_text(size = 12,face="bold"),
                                      axis.title.x = element_blank(),
                                      axis.text = element_text(size = 10,face="bold"))+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(xx$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(xx$niche)*1.1, x = -Inf+1, xend = -Inf+1)+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(xx$group))+0.7/2,
             y = Inf, yend = Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(xx$niche)*1.1, x = Inf+1, xend = Inf+1)

  if (panel.border=="axis") p<-p+theme(aspect.ratio = 1,
                                       title =  element_text(size = 7,face="bold"),
                                       axis.title = element_text(size = 12,face="bold"),
                                       axis.title.x = element_blank(),
                                       axis.text = element_text(size = 10,face="bold"))+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(xx$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(xx$niche)*1.1, x = -Inf+1, xend = -Inf+1)

  if (panel.border=="normal") p<-p+ylim(0,max(xx$niche)*1.1)+
    theme(
    aspect.ratio = 1,
    title =  element_text(size = 7,face="bold"),
    axis.title = element_text(size = 12,face="bold"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10,face="bold"),
    panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5)
  )
  if (panel.border=="none") p<-p

  if (axis.x.angle.adjust) p<-p+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
  if (!is.null(mytheme)) p<-p+mytheme

  if (!is.null(export_path)) {
    file2=paste(export_path,"/Ecological niche width/","Allgroup_niche breadth.pdf",sep = "")

    ggsave(file2,
           units = "cm",
           width = 21/3,
           height = 21*p$theme$aspect.ratio/3,
           p)
    cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment has been exported to ","/",export_path,"",sep = "","\n")
  } else {
    cat("\n","Relationship between mean abundance of OTUs and ecological niche width in ",names(spec_gen)[tt]," treatment was not exported to local path","",sep = "","\n")
  }

  return(p)
}

"fun.outlier" <- function(x,time.iqr=1.5) {
  outlier.low <- quantile(x,probs=c(0.25))-IQR(x)*1.5
  outlier.high <- quantile(x,probs=c(0.75))+IQR(x)*1.5
  return(list(outlier.low=outlier.low,outlier.high=outlier.high))
}

"net_visulal_chord" <- function(df,grid.color,order.sector) {

  df[,c("phylum","genus")]<-stringr::str_split_fixed(df$Var1,pattern=";",2)
  df.plot<-subset(df, select = c(3,4,2))

  circlize::circos.clear()
  par(family="serif",mar=c(2,2,2,2))
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
                          niceFacing = FALSE, adj = c(0, 0), cex = 0.5)
  }, bg.border = NA)

  p1<-recordPlot()
  circlize::circos.clear()

  message("\n","The chord plot need to be saved manually.")

  return(p1)
}


"plotMicrochatAssemblySourceChord"  <- function(spec_gen,
                                                specified.thres=NULL,
                                                nrow=2,
                                                color_phylum,
                                                color_genus) {

  sepc.info<-spec_gen$spec_gen
  message("specified.thres could be a vector, whose amount is equal with treatment.")
  if (is.null(specified.thres)) high.thres<-spec_gen$high.thres else high.thres<-specified.thres
  sepc.info1<-list()
  sepc.info2<-list()
  for (tt in 1:length(sepc.info)) {
    sepc.info1[[tt]]<-sepc.info[[tt]][which(sepc.info[[tt]]$niche>=high.thres[tt]),]
    sepc.info2[[tt]]<-sepc.info[[tt]][which(sepc.info[[tt]]$niche<1.5),]
  }

  sepc.info<-c(sepc.info1,sepc.info2)
  dd.data<-data.frame()
  for (tk in 1:length(sepc.info)) {
    gma<-names(sepc.info)[tk]
    dd.dat<-sepc.info[[tk]]
    dd.dat$group<-gma
    dd.data<-rbind(dd.data,dd.dat)
  }
  phy.num<-length(unique(dd.data$Phylum))
  gen.num<-length(unique(dd.data$Genus))

  while (length(color_phylum)<phy.num) {
    color_phylum<-c(color_phylum,color_phylum)
  }
  color_phylum<-color_phylum[1:phy.num]

  while (length(color_genus)<gen.num) {
    color_genus<-c(color_genus,color_genus)
  }
  color_genus<-color_genus[1:gen.num]

  names(color_phylum)<-unique(dd.data$Phylum)
  names(color_genus)<-unique(dd.data$Genus)
  color.comb<-c(color_phylum,color_genus)

  if (!is.null(nrow)) ncol=nrow else ncol=ncol_layout(length(sepc.info))%>%length()
  par(mfrow=c(ncol,length(sepc.info)/ncol))
  for (tt in 1:length(sepc.info)) {
    sepc.info[[tt]]$com<-paste(sepc.info[[tt]]$Phylum,sepc.info[[tt]]$Genus,sep = ";")
    df<-table(sepc.info[[tt]]$com)%>%data.frame()
    df[,c("phylum","genus")]<-stringr::str_split_fixed(df$Var1,pattern=";",2)
    df.plot<-subset(df, select = c(3,4,2))

    order.sector<- c(unique(df.plot$phylum)[length(unique(df.plot$phylum)):1],
                     unique(df.plot$genus))

    grid.color<-color.comb[match(order.sector,names(color.comb))]
    net_visulal_chord(df,grid.color,order.sector)
  }
  pp<-recordPlot()
  message("Generalists from all treatments was listed in the first row, followed by the specialists listed in the second row!!!")
  return(pp)
}



"plotMicrochatAssemblySourceCircle" <- function(spec_gen,
                                                edge.curved=0.5,
                                                nrow=2,
                                                sel.layout,
                                                specified.thres=NULL,
                                                color_phylum,
                                                color_genus) {
  sepc.info<-spec_gen$spec_gen
  message("specified.thres could be a vector, whose amount is equal to treatment.")
  if (is.null(specified.thres)) high.thres<-spec_gen$high.thres else high.thres<-specified.thres
  sepc.info1<-list()
  sepc.info2<-list()
  for (tt in 1:length(sepc.info)) {
    sepc.info1[[tt]]<-sepc.info[[tt]][which(sepc.info[[tt]]$niche>=high.thres[tt]),]
    sepc.info2[[tt]]<-sepc.info[[tt]][which(sepc.info[[tt]]$niche<1.5),]
  }

  sepc.info<-c(sepc.info1,sepc.info2)

  dd.data<-data.frame()
  for (tk in 1:length(sepc.info)) {
    gma<-names(sepc.info)[tk]
    dd.dat<-sepc.info[[tk]]
    dd.dat$group<-gma
    dd.data<-rbind(dd.data,dd.dat)
  }
  phy.num<-length(unique(dd.data$Phylum))
  gen.num<-length(unique(dd.data$Genus))

  while (length(color_phylum)<phy.num) {
    color_phylum<-c(color_phylum,color_phylum)
  }
  color_phylum<-color_phylum[1:phy.num]

  while (length(color_genus)<gen.num) {
    color_genus<-c(color_genus,color_genus)
  }
  color_genus<-color_genus[1:gen.num]

  names(color_phylum)<-unique(dd.data$Phylum)
  names(color_genus)<-unique(dd.data$Genus)
  color.comb<-c(color_phylum,color_genus)


  roles<-rep(c('Generalist',"Specialist"),each=length(sepc.info1))
  if (!is.null(nrow)) ncol=nrow else ncol=ncol_layout(length(sepc.info))%>%length()
  par(mfrow=c(ncol,length(sepc.info)/ncol),mar=c(5,5,5,5),family="serif",font=4)

  for (tt in 1:length(sepc.info)) {
    gname<-c(names(spec_gen$spec_gen),names(spec_gen$spec_gen))[tt]
    sel.role<-roles[tt]
    sepc.info[[tt]]$com<-paste(sepc.info[[tt]]$Phylum,sepc.info[[tt]]$Genus,sep = ";")
    df<-table(sepc.info[[tt]]$com)%>%data.frame()
    df[,c("phylum","genus")]<-stringr::str_split_fixed(df$Var1,pattern=";",2)
    df.plot<-subset(df, select = c(3,4,2))

    ddf.plot<-subset(df.plot,select = c(phylum,Freq))%>%group_by(phylum)%>%summarise_all(sum)
    ddf.plot<-subset(ddf.plot,select = c(1,1,2))
    colnames(ddf.plot)[2]<-"genus"
    df.plot<-rbind(ddf.plot,df.plot)
    df.plot<-df.plot[order(df.plot$phylum),]
    df.mat<-subset(df.plot,select=c(phylum,genus))%>%as.matrix()

    gg<-igraph::graph_from_edgelist(df.mat[,1:2],directed = FALSE)

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
    edge.start <- igraph::ends(gg, es = igraph::E(gg), names = FALSE)
    if (sel.layout=="circle") loop.angle <- ifelse(coords_scale[igraph::V(gg), 1] > 0,
                                                   -atan(coords_scale[igraph::V(gg), 2]/coords_scale[igraph::V(gg),
                                                                                                     1]), pi - atan(coords_scale[igraph::V(gg), 2]/coords_scale[igraph::V(gg),
                                                                                                                                                                1]))

    if (sel.layout=="circle") igraph::E(gg)$loop.angle[which(edge.start[, 2] == edge.start[,
                                                                                           1])] <- loop.angle[edge.start[which(edge.start[,
                                                                                                                                          2] == edge.start[, 1]), 1]]


    E(gg)$weight<-df.plot$Freq
    E(gg)$width <- E(gg)$weight
    igraph::V(gg)$color <- color.comb[match(igraph::V(gg)$name,names(color.comb))]
    igraph::V(gg)$frame.color <- color.comb[match(igraph::V(gg)$name,names(color.comb))]
    igraph::E(gg)$color <- grDevices::adjustcolor(igraph::V(gg)$color[edge.start[,1]], 0.6)

    "radian.rescale" <- function(x, start = 0, direction = 1) {
      c.rotate <- function(x) (x + start)%%(2 * pi) * direction
      c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }
    if (sel.layout=="circle") label.locs <- radian.rescale(x = 1:length(igraph::V(gg)),
                                                           direction = -1, start = 0)

    vertex.size.sel<-E(gg)$weight
    remove.thres<-mean(vertex.size.sel)%>%as.integer()
    V(gg)$name[substr(V(gg)$name,start = 1,stop = 3)=="g__" & E(gg)$weight<=remove.thres]<-""
    vertex.size.sel[vertex.size.sel<3]<-2.9
    vertex.size.sel[vertex.size.sel<6 & vertex.size.sel>=3]<-5.9
    vertex.size.sel[vertex.size.sel<9 & vertex.size.sel>=6]<-8.9
    vertex.size.sel[vertex.size.sel<12 & vertex.size.sel>=9]<-11.9
    vertex.size.sel[vertex.size.sel>=12]<-12
    V(gg)$size<-vertex.size.sel
    E(gg)$width<-vertex.size.sel
    message("Node size range: 2.9--5.9--8.9--11.9--12")
    if (sel.layout=="circle") plot(gg,edge.curved=edge.curved,
                                   rescale = TRUE,
                                   layout=coords,
                                   vertex.size=vertex.size.sel,
                                   vertex.shape="circle",
                                   vertex.label.angle=loop.angle,
                                   vertex.label.degree = label.locs,
                                   vertex.label.dist = 4,
                                   vertex.label=vertex_attr(gg)$name,
                                   vertex.frame.color = "white",
                                   vertex.label.family = "serif",
                                   edge.label.family = "serif")

    if (sel.layout!="circle") plot(gg,edge.curved=edge.curved,
                                  rescale = TRUE,
                                  layout=coords,
                                  vertex.size=vertex.size.sel,
                                  vertex.shape="circle",
                                  vertex.label.dist = 3,
                                  vertex.label=vertex_attr(gg)$name,
                                  vertex.frame.color = "white",
                                  vertex.label.family = "serif",
                                  edge.label.family = "serif")
    title(main = paste0(gname,"--",sel.role))

  message("Genus label that appears once has been hidden.")
  }
  pp<-recordPlot()
  return(pp)
}

"calc_NCM" <- function(spp) {

  spp<-t(spp)
  N <- mean(apply(spp, 1, sum))
  p.m <- apply(spp, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  d = 1/N
  m.fit <- minpack.lm::nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
  m.fit
  m.ci <- confint(m.fit, 'm', level=0.95)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
  pred.ci <- Hmisc::binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  Rsqr
  return(list(p=p,freq=freq,freq.pred=freq.pred,pred.ci=pred.ci,Rsqr=Rsqr,m.fit=m.fit))
}


"calcMicrochatAssemblyNCM" <- function(submchat,
                                       export_path="cs2/microbial assembly analysis") {
  export_path<-paste(export_path,"/data_microbiome/microbial assembly analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,"/Neutral community model",sep = ""), recursive = TRUE)

  if (class(submchat)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  sample.size=submchat$otu_table%>%colSums()
  if (length(unique(sample.size))==1) sample.size<-sample.size[1] else sample.size<-max(sample.size)
  otutab<-submchat$otu_table
  split_otux<-calc_indivGroup_abun(otutab)

  split_otux.mean<-(sapply(split_otux, rowMeans)/sample.size*100)%>%data.frame()
  split_otux.mean$Whole<-rowMeans(split_otux.mean)
  split_otux.mean<-split_otux.mean%>%tibble::rownames_to_column(var = "name")

  plot.data<-lapply(split_otux, calc_NCM)
  plot.data.all<-calc_NCM(otutab)
  plot.data[["Whole"]]<-plot.data.all

  for (kt in 1:length(plot.data)) {
    sel.abun<-split_otux.mean[,c(1,kt+1)]
    colnames(sel.abun)[2]<-"abun"
    plot.data[[kt]]$abundance<-sel.abun
    plot.data[[kt]]$sample.size<-sample.size
  }


  for (tt in 1:length(plot.data)) {
    gname<-names(plot.data)[tt]
    filex<-paste(export_path,"/Neutral community model/",gname,sep = "")
    dir.create(filex, recursive = TRUE)
    plot.dat<-plot.data[[tt]]
    p<-plot.dat$p
    freq<-plot.dat$freq
    freq.pred<-plot.dat$freq.pred

    if (!is.null(export_path)) {
      write.csv(p, file = paste(filex,"/p.csv",sep = ""))
      write.csv(freq, file = paste(filex,"/freq.csv",sep = ""))
      write.csv(freq.pred, file = paste(filex,"/freq.pred.csv",sep = ""))
      cat("\n","Neutral community model in ",gname,"  has been exported to ","/",filex,"",sep = "","\n")
    } else {
      cat("\n","Neutral community model in ",gname,"  was not exported to local path","",sep = "","\n")
    }
  }


  message("Fitting model parameters with non-linear least squares (NLS).")
  message("Notes:")
  message("  p: mean relative abundance")
  message("  freq: the observed value of occurrence frequency")
  message("  freq.pred: the predicted value of occurrence frequency--fitted value of NCM")

  return(plot.data)
}

"plotMicrochatAssemblyNCM" <- function(plot.ddat,
                                       ncol=3,
                                       xlabname=NULL,
                                       nrow_layout=NULL,
                                       point.size=1,
                                       add.border=FALSE,
                                       color_background=NA,
                                       color_pie.text="white",
                                       color_NCM=c('black','#A52A2A','#29A6A6'),
                                       color_line=c("blue","blue","blue"),
                                       mytheme=NULL,
                                       export_path="cs2/microbial assembly analysis") {
  export_path<-paste(export_path,"/data_microbiome/microbial assembly analysis",sep = "")

  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,"/Neutral community model",sep = ""), recursive = TRUE)

  xmax<-sapply(plot.ddat, function(x){
    max(log10(x$p))
  })%>%max()
  xmin<-sapply(plot.ddat, function(x){
    min(log10(x$p))
  })%>%min()

if (!is.null(xlabname)) names(plot.ddat)[1:(length(plot.ddat)-1)]<-xlabname

  pp<-list()
  for (tk in 1:length(plot.ddat)) {

    pE<-plot.ddat[[tk]]$p
    freq<-plot.ddat[[tk]]$freq
    freq.pred<-plot.ddat[[tk]]$freq.pred
    pred.ci<-plot.ddat[[tk]]$pred.ci
    m.fit<-plot.ddat[[tk]]$m.fit
    Rsqr<-plot.ddat[[tk]]$Rsqr
    abun<-plot.ddat[[tk]]$abundance
    N<-plot.ddat[[tk]]$sample.size

    bacnlsALL <-data.frame(pE,freq,freq.pred,pred.ci[,2:3])
    bacnlsALL$inter.col<-"NCM"
    bacnlsALL$inter.col[which(bacnlsALL$freq <= bacnlsALL$Lower)]<-'Low'
    bacnlsALL$inter.col[which(bacnlsALL$freq >= bacnlsALL$Upper)]<-'High'
    col<-color_NCM
    names(col)<-c("NCM","Low","High")
    bacnlsALL$pE<-log10(bacnlsALL$pE)

    bacnlsALL<-bacnlsALL%>%tibble::rownames_to_column(var = "name")
    bacnlsALL<-merge(bacnlsALL,abun,by="name")
    #length(unique(bacnlsALL$inter.col))
    bacnlsALL.x<-bacnlsALL
    bacnlsALL.x$inter.col<-factor(bacnlsALL.x$inter.col,levels = c("Low","NCM","High"))
    abun.all<-subset(bacnlsALL.x,select = c(inter.col,abun))%>%
      group_by(inter.col)%>%
      summarise_all(sum)
    abun.all<-abun.all[order(abun.all$abun,decreasing = TRUE),]

    barph<-ggplot(abun.all,aes(x="",y=abun)) +
      geom_bar(stat = "identity",show.legend = FALSE,
               aes(fill=inter.col),width = 0.42)+
      geom_text(aes(label=ifelse(abun>3.9,paste(round(abun,0),"%",sep = ""),"")),
                color=color_pie.text,family="serif",size=2,
                position = position_stack(vjust = 0.5))+
      scale_fill_manual(values = col)+
      coord_flip()+
      theme(text = element_text(family = "serif"),
            plot.background = element_rect(fill = NA,color = NA),
            legend.position = "none")+
      theme_void()

    pp[[tk]]<-ggplot(bacnlsALL,aes(x=pE))+
      geom_point(aes(y=freq,color=inter.col),size=point.size)+
      geom_line(aes(y=freq.pred),
                color=color_line[1])+
      geom_line(aes(y=Lower),
                color=color_line[2],
                linetype=2)+
      geom_line(aes(y=Upper),
                color=color_line[3],
                linetype=2)+
      scale_color_manual(values = col,name="NCM fitting")+
      labs(x='Mean Relative Abundance (log10)',
           y='Frequency of Occurance',title = names(plot.ddat)[tk])+
      annotate(geom = "text",
               x = xmin*0.85, y = 0.85,
               family="serif",size=2.5,
               label = paste0("R2=",round(Rsqr,3),"\n",
                             "Nm=",round(coef(m.fit)*N),"\n",
                             "N=",N,"\n","m=",round(coef(m.fit),3)))+
      theme(text = element_text(family = "serif"),
            legend.position = "none",
            #plot.margin = margin(r = 20, l = 20, unit = "pt"),
            panel.background = element_rect(fill = color_background),
            aspect.ratio = 1)+
      scale_x_continuous(limits = c(xmin,xmax))

    pp[[tk]]<-pp[[tk]]+theme(aspect.ratio = 1,
                 title =  element_text(size = 8,face="bold"),
                 axis.title = element_text(size =10,face="bold"),
                 axis.text = element_text(size = 8,face="bold"),
                 plot.background = element_rect(fill = NA,color = NA),
                 panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5))
    if (!is.null(mytheme)) pp[[tk]]<-pp[[tk]]+mytheme

    if (add.border) pp[[tk]]<- pp[[tk]]+ theme(panel.border = element_rect(fill=NA))

    node_phylum<-table(bacnlsALL$inter.col)%>%data.frame()
    node_phylumx<-node_phylum
    node_phylumx<-node_phylumx%>%
      group_by(Var1)%>%
      summarise_all(sum)
    node_phylumx<-node_phylumx[order(node_phylumx$Freq,decreasing = TRUE),]
    node_phylumx$Var1<-factor(node_phylumx$Var1,levels =rev(unique(node_phylumx$Var1)) )

    barpv<-ggplot(node_phylumx,aes(x="",y=Freq)) +
      geom_bar(stat = "identity",show.legend = FALSE,
               aes(fill=Var1),width = 0.42)+
      geom_text(aes(label=ifelse(Freq/sum(Freq)*100>3.9,
                                 paste(round(Freq/sum(Freq)*100,0),"%",sep = ""),
                                 "")),
                color=color_pie.text,family="serif",angle=270,size=2,
                position = position_stack(vjust = 0.5))+
      scale_fill_manual(values = col)+
      theme(text = element_text(family = "serif"),
            plot.background = element_rect(fill = NA,color = NA),
            legend.position = "none")+
      theme_void()

    piep<-ggplot(node_phylum, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', width = 0.5) +
      coord_polar(theta = "y") +
      scale_fill_manual(limits = names(col), values = col) +

      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.05,
                    label=ifelse(Freq/sum(node_phylum$Freq)*100>3.9,
                       paste(round(Freq/sum(node_phylum$Freq)*100,0),"%",sep = ""),
                       "")),
                size=2,colour=color_pie.text,family = "serif")+
      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.35,label=Var1),fontface="bold",
                size=2,color="black",family = "serif")+
      theme(panel.grid = element_blank(),
            plot.background = element_rect(fill = NA,color = NA),
            axis.text.x = element_blank()) +
      labs(x = '', y = '')+
        theme_void()+theme(legend.position = "none")

    (xxmin<-(xmin-xmax)/2)
    pp[[tk]]<- pp[[tk]] + annotation_custom(grob=ggplotGrob(piep),
                                 ymin = -0.2, ymax = 0.6,
                                 xmin=xxmin-0.3,
                                 xmax=xmax+0.3)

    ### add high layout
    #pp[[tk]]<- pp[[tk]] + annotation_custom(grob=ggplotGrob(barpv),ymin = 1, ymax = 1.15,xmin=xmin-0.525,xmax=xmax+0.525)

    ### add right layout
    pp[[tk]]<- pp[[tk]] + annotation_custom(grob=ggplotGrob(barpv),
                                            ymin = -0.105, ymax = 1.105,
                                            xmin=xmax,
                                            xmax=xmax+0.75)

    pp[[tk]]<- pp[[tk]] + annotation_custom(grob=ggplotGrob(barph),
                                            ymin = -0.1, ymax = 0.05,
                                            xmin=xmin-0.525,
                                            xmax=xmax+0.525)




    if (!is.null(export_path)) {
      filex<-paste(export_path,"/Neutral community model/",names(plot.ddat)[tk],"_NCM.pdf",sep = "")
      ggsave(filex,
             units = "cm",
             width = 7,
             height = 7*pp[[tk]]$theme$aspect.ratio,
             pp[[tk]])

      cat("\n","Neutral community model in ",names(plot.ddat)[tk],"  has been exported to ","/",filex,"",sep = "","\n")
    } else {
      cat("\n","Neutral community model in ",names(plot.ddat)[tk],"  was not exported to local path","",sep = "","\n")
    }

  }

  if (!is.null(ncol)){
    width.sel<-ncol
    if ((length(pp)/width.sel)%%1==0) {
      height.sel<-(length(pp)/width.sel)%>%as.integer()
    } else {
      height.sel<-(length(pp)/width.sel)%>%as.integer()+1
    }
  } else {
    width.sel<-ncol_layout(length(pp))[[1]]%>%length()
    height.sel<-ncol_layout(length(pp))%>%length()
  }
  if (!is.null(nrow_layout)) xxxxpx<-nrow_layout
  mp<-patchwork::wrap_plots(pp,ncol = width.sel)
  if (!is.null(export_path)) {
    filex<-paste(export_path,"/Neutral community model/","Muti-NCM.pdf",sep = "")
    rownum<-length(pp)/ncol
    if (rownum>as.integer(rownum)) nrows<-as.integer(rownum)+1 else nrows<-as.integer(rownum)
    ggsave(filex,
           units = "cm",
           width = 7*ncol,height = 7*nrows,
           mp)

    cat("\n","Neutral community model in whole data"," has been exported to ","/",filex,"",sep = "","\n")
  } else {
    cat("\n","Neutral community model "," was not exported to local path","",sep = "","\n")
  }

  message("R2 represents the overall fit of the neutral community model. A higher R2 indicates a closer approximation to the neutral model, meaning that the construction of the community is more influenced by stochastic processes and less influenced by deterministic processes")
  message("R2代表了中性群落模型的整体拟合优度，R2越高表明越接近中性模型，即群落的构建受随机性过程的影响越大，受确定性过程的影响越小")
  message("N描述了每个样本中OTU的总丰度；m量化了群落层面的迁移率（migration rate），m值越小说明整个群落中物种扩散越受限制，反之m值越高则表明物种受到扩散限制越低；Nm是群落规模（N）与迁移率（m）的乘积，量化了对群落之间扩散的估计，决定了发生频率和区域相对丰度之间的相关性。")
  message("生态位理论（niche-based theories）和中性理论（neutral-based theories）构成了理解微生物群落组装的两个重要且互补的机制。生态位过程认为微生物群落是由确定性的非生物因素（环境因素如pH、温度等）和生物因素（物种相互作用如竞争、捕食等）所决定，归因于微生物不同的栖息地偏好和适应能力。中性理论则认为，随机过程，如出生、死亡、迁移、物种形成和有限扩散，塑造了微生物群落结构，并假设微生物在类群的损失和收益之间呈现一种随机平衡。")
  message("Neutral model applied to assess the effects of random dispersal and ecological drift on the assembly of microbial community.")
  message("m indicates the estimated migration rate. The solid blue lines indicate the best fit to the neutral model and dashed blue lines represent 95% confidence intervals around the model prediction.")
  return(mp)
}



"calBNTIx" <- function(data = match.phylo.otu,thres=100,nworker=4){
  require(picante)
  require(tidyverse)
  message("6 samples with 999 (thres) will cost about 100 seconds. 100 (thres) was recommended.")
  set.seed(123)
  beta.mntd.weighted = as.matrix(comdistnt(t(data$data),
                                           cophenetic(data$phy),
                                           abundance.weighted=T))

  rand.weighted.bMNTD.comp = array(c(-thres),dim=c(ncol(beta.mntd.weighted),
                                                   ncol(beta.mntd.weighted),thres))

  dim(rand.weighted.bMNTD.comp)

  "rd.w.mntd" <- function(rep,rand.weighted.bMNTD.comp,data) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(picante::comdistnt(t(data$data),
                                                                   picante::taxaShuffle(cophenetic(data$phy)),
                                                                   abundance.weighted=T,exclude.conspecifics = F))

    return(rand.weighted.bMNTD.comp)
  }

  if (nworker!=1) {
    c1 <- try(parallel::makeCluster(nworker, type = "PSOCK"))
    if (class(c1)[1] == "try-error") {
      c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))
    }
    if (class(c1)[1] == "try-error") {
      c1 <- try(parallel::makeCluster(nworker, setup_strategy = "sequential"))
    }
    message("Now randomizing by parallel computing. Begin at ",
            date(), ". Please wait...")
    rand.weighted.bMNTD.comp = parallel::parLapply(c1, 1:thres, rd.w.mntd,
                                                   rand.weighted.bMNTD.comp,
                                                   data)
    parallel::stopCluster(c1)
  }

  rand.weighted.bMNTD.compx = array(c(-thres),dim=c(ncol(beta.mntd.weighted),
                                                    ncol(beta.mntd.weighted),thres))
  for (i in 1:length(rand.weighted.bMNTD.comp)) {
    rand.weighted.bMNTD.compx[,,i]<-rand.weighted.bMNTD.comp[[i]][,,i]
  }
  rand.weighted.bMNTD.comp<-rand.weighted.bMNTD.compx


  weighted.bNTI <- matrix(c(NA),nrow=ncol(beta.mntd.weighted),
                          ncol=ncol(beta.mntd.weighted))

  for (columns in 1:(ncol(beta.mntd.weighted)-1)) {
    for (rows in (columns+1):ncol(beta.mntd.weighted)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,]
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
      rm("rand.vals")

    }
  }
  rownames(weighted.bNTI) <- colnames(data$data)
  colnames(weighted.bNTI) <- colnames(data$data)

  betaNTI <- weighted.bNTI %>%
    as_tibble(rownames = "sample") %>%
    reshape2::melt(id.vars = "sample") %>%
    dplyr::filter(value!="NA") %>%
    mutate(process = if_else(value > 2, "Deter_HeS",
                             if_else(value >= -2,"Stochastic","Deter_HoS")))

  return(betaNTI)
}

"calcMicrochatAssemblyBNT" <- function(submchat,
                                       bootstrap.thres=100,
                                       nworker=6) {
  test_tree<-submchat$tree
  test_otu<-submchat$otu_table%>%t()%>%data.frame()
  otutab<-calc_indivGroup_abun(submchat$otu_table)

  bnt<-data.frame()
  for (ttk in 1:length(otutab)) {
    gname<-names(otutab)[ttk]
    test_otu<-otutab[[ttk]]%>%t()%>%data.frame()
    test_otu<-test_otu[,-which(colSums(test_otu)==0)]
    prune_tree <- picante::prune.sample(test_otu,test_tree)#select tree data according to otu data
    phydist<-cophenetic(prune_tree)

    t_test_otu <- t(test_otu) %>%
      data.frame()
    match_phylo_otu <- match.phylo.data(prune_tree,t_test_otu)

    a<-Sys.time()
    βNTI <- calBNTIx(match_phylo_otu,thres=bootstrap.thres,nworker=nworker)
    βNTI
    b<-Sys.time()
    costt<-(b-a)
    class(costt)
    message(gname,": βNTI calculation cost ", b-a,".")
    βNTI$group<-gname

    {
      bnt<-rbind(bnt,βNTI)
    }
  }
  return(bnt)
}

"plotMicrochatAssemblyBNT" <- function(BNT,
                                       geom=c("barplot","pie","chord"),
                                       pie.layout.nrow=1,
                                       color_group=colorCustom(4,pal = "gygn"),
                                       color_model=colorCustom(2,pal = "gygn"),
                                       xlabname=c("CS0","CS100","CA100","CD100"),
                                       export_path="cs2/microbial assembly analysis") {

  geom<-match.arg(geom)
  suppressMessages(library(circlize))
  export_path<-paste(export_path,"/data_microbiome/microbial assembly analysis",sep = "")

  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,"/Beta nearest taxon index (βNTI)",sep = ""), recursive = TRUE)

  BNT %>%
    dplyr::group_by(group) %>%
    dplyr::count(process) %>%
    dplyr::mutate(prop=round(n/sum(n)*100,2)) %>%
    data.frame()->dat

  explan.index<-data.frame(x="*βNTI>2 indicates 超过预期的系统发育周转，解释为确定性过程的异质选择(variable selection of deterministic process)",
                           y="*βNTI<-2 indicates 低于预期的系统发育周转，解释为确定性过程的同质选择(homogenuous selection of deterministic process)",
                           z="*-2<=βNTI<=2 indicates 观察到的系统发育组成差异是随机过程的结果")
  explan.index<-cbind(explan.index,
                      data.frame(matrix(NA,
                                        nrow=1,
                                        ncol=ncol(BNT)-ncol(explan.index))))
  colnames(explan.index)<-colnames(BNT)
  BNTx<-rbind(BNT,explan.index)
  filex<-paste(export_path,"/Beta nearest taxon index (βNTI)/",sep = "")
  if (!is.null(export_path)) {
    write.csv(BNTx, file = paste(filex,"Beta nearest taxon index (βNTI).csv",sep = ""))
    write.csv(dat, file = paste(filex,"βNTI proportion.csv",sep = ""))
    cat("\n","Beta nearest taxon index (βNTI)"," has been exported to ","/",filex,"",sep = "","\n")
  } else {
    cat("\n","Beta nearest taxon index (βNTI)"," was not exported to local path","",sep = "","\n")
  }

  dat$group<-factor(dat$group,levels = unique(BNT$group))
  dat$process<-factor(dat$process,levels = unique(BNT$process))
  dat<-dat[order(dat$group,dat$process),]

  names(color_model)<-unique(dat$process)

  ###barplot
  if (geom=="barplot") {
    bntp<-ggplot(data=dat,aes(group,prop,fill=process))+
      geom_col(position = "stack",
               width = 0.7,linetype="solid")+
      geom_text(aes(label=paste(prop,"%",sep = "")),family="serif",
                position = position_stack(vjust = 0.5))+
      labs(y="Model dominant proportion (%)")+
      scale_fill_manual(values = color_model,name="Model")+
      scale_x_discrete(labels=xlabname)+
      theme(aspect.ratio = 1,
            axis.title.x = element_blank(),
            legend.position = "bottom",
            text = element_text(family = "serif"))
  }

  ###pie

  if (geom=="pie") {
    pp<-list()
    for (gp in unique(dat$group)) {
      datx<-dat[which(dat==gp),]
      gorder<-match(gp,unique(dat$group))
      p<-ggplot(data=datx,aes(x="",prop,fill=process))+
        geom_bar(stat = "identity",
                 width = 0.4)+
        coord_polar(theta = "y") +
        scale_fill_manual(values = color_model,name="Model") +
        geom_text(aes(y=rev(rev(prop/2)+c(0,cumsum(rev(prop))[-length(prop)])),
                      x= 1.05,label=scales::percent(prop/sum(prop))),
                  size=4,colour="white",family = "serif")+
        geom_text(aes(y=rev(rev(prop/2)+c(0,cumsum(rev(prop))[-length(prop)])),
                      x= 1.35,label=process),
                  size=4,color="black",family = "serif")+
        labs(x = '', y = '',title = xlabname[gorder])+
        theme_void()+
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              plot.title  = element_text(hjust = 0.5,family = "serif"))


      pp[[gp]]<-p
    }
    bntp<-patchwork::wrap_plots(pp,nrow = pie.layout.nrow)
  }

  if (geom=="chord") {
    ###chord
    df.plot<-dat[,c(1:2,4)]
    order.sector<- c(unique(df.plot$process)[length(unique(df.plot$process)):1],
                     unique(df.plot$group))

    color_phylum<-color_model
    color_genus<-color_group
    names(color_phylum)<-unique(df.plot$process)
    names(color_genus)<-unique(df.plot$group)
    color.comb<-c(color_phylum,color_genus)
    grid.color<-color.comb[match(order.sector,names(color.comb))]

    circlize::circos.clear()
    par(family="serif",mar=c(2,2,2,2))
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
                            facing = "bending.inside",family = par("serif"),
                            niceFacing = FALSE, adj = c(0.5, -0.5), cex = 1)
    }, bg.border = NA)
    suppressWarnings(warning("circos.track"))
    bntp<-recordPlot()
  }

  if (!is.null(export_path)) {
    filex<-paste(export_path,"/Beta nearest taxon index (βNTI)/βNTI proportion_",geom,".pdf",sep = "")
    if (geom=="chord") {
     message("Please save the plot mannully!!!")
    } else  ggsave(filex, bntp)
    cat("\n","Beta nearest taxon index (βNTI) proportion","  has been exported to ","/",filex,"",sep = "","\n")
  } else {
    cat("\n","Beta nearest taxon index (βNTI) proportion","  was not exported to local path","",sep = "","\n")
  }
  return(bntp)
}


"calcMicrochatAssemblyNST" <- function(submchat,
                                       nworker=4,
                                       dist.method="bray") {

  test_otu<-submchat$otu_table%>%t()%>%data.frame()
  group<-group_generate(submchat$otu_table)%>%tibble::column_to_rownames(var = "sample")

  tNST_test <- NST::tNST(comm=test_otu, group=group,
                         dist.method=dist.method, abundance.weighted=TRUE, rand=999,
                         output.rand=TRUE, nworker=nworker, LB=FALSE, null.model="PF",
                         between.group=FALSE, SES=TRUE, RC=TRUE)


  tNST_test[[1]]#indlude NST(i.e. ST.ij.bray) and RC(i.e. RC.bray)
  tNST_test[[2]]#tNST (i.e., NST.i.bray)
  tNST_test[[3]] %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(value = mean(NST.ij.bray))


  nst<-tNST_test
  return(nst)
}




"plotMicrochatAssemblyNST" <- function(nst,
                                       geom=c("barplot","pie","chord"),
                                       pie.layout.nrow=1,
                                       color_group=colorCustom(4,pal = "gygn"),
                                       color_model=colorCustom(2,pal = "gygn"),
                                       xlabname=c("CS0","CS100","CA100","CD100"),
                                       export_path="cs2/microbial assembly analysis") {

  geom<-match.arg(geom)
  suppressMessages(library(circlize))
  export_path<-paste(export_path,"/data_microbiome/microbial assembly analysis",sep = "")

  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,"/Normalized stochasticity ratio (NST)",sep = ""), recursive = TRUE)


  nst[[2]] %>%
    subset(select=c(1,4)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(deter=1-NST.i.bray) %>%
    data.frame()->dat
  colnames(dat)[2]<-"Stochastic"
  colnames(dat)[3]<-"Deterministic"

  explan.index<-data.frame(x="*NST>0.5 解释为随机性过程(Stochastic process)",
                           y="*NST<0.5 解释为确定性过程(Deterministic process)")
  explan.index<-cbind(explan.index,
                      data.frame(matrix(NA,
                                        nrow=1,
                                        ncol=ncol(nst[[2]])-ncol(explan.index))))
  colnames(explan.index)<-colnames(nst[[2]])
  nst[[2]]<-rbind(nst[[2]],explan.index)
  filex<-paste(export_path,"/Normalized stochasticity ratio (NST)/",sep = "")
  if (!is.null(export_path)) {
    write.csv(nst[[2]], file = paste(filex,"Normalized stochasticity ratio (NST).csv",sep = ""))
    write.csv(dat, file = paste(filex,"NST proportion.csv",sep = ""))
    write.csv(nst[[3]], file = paste(filex,"Normalized stochasticity ratio (NST)_ori.csv",sep = ""))
    cat("\n","Normalized stochasticity ratio (NST)"," has been exported to ","/",filex,"",sep = "","\n")
  } else {
    cat("\n","Normalized stochasticity ratio (NST)"," was not exported to local path","",sep = "","\n")
  }

  dat<-reshape2::melt(dat)
  colnames(dat)[2]<-"process"
  dat$prop<-sprintf("%0.2f",round(dat$value*100,3))%>%as.numeric()

  dat$group<-factor(dat$group,levels = unique(nst[[2]]$group))
  dat$process<-factor(dat$process,levels = unique(dat$process))
  dat<-dat[order(dat$group,dat$process),]

  names(color_model)<-unique(dat$process)

  ###barplot
  if (geom=="barplot") {
    bntp<-ggplot(data=dat,aes(group,prop,fill=process))+
      geom_col(position = "stack",
               width = 0.7,linetype="solid")+
      geom_text(aes(label=paste(prop,"%",sep = "")),family="serif",
                position = position_stack(vjust = 0.5))+
      labs(y="Model dominant proportion (%)")+
      scale_fill_manual(values = color_model,name="Model")+
      scale_x_discrete(labels=xlabname)+
      theme(aspect.ratio = 1,
            axis.title.x = element_blank(),
            legend.position = "bottom",
            text = element_text(family = "serif"))
  }

  ###pie

  if (geom=="pie") {
    pp<-list()
    for (gp in unique(dat$group)) {
      datx<-dat[which(dat==gp),]
      gorder<-match(gp,unique(dat$group))
      p<-ggplot(data=datx,aes(x="",prop,fill=process))+
        geom_bar(stat = "identity",
                 width = 0.4)+
        coord_polar(theta = "y") +
        scale_fill_manual(values = color_model,name="Model") +
        geom_text(aes(y=rev(rev(prop/2)+c(0,cumsum(rev(prop))[-length(prop)])),
                      x= 1.05,label=scales::percent(prop/sum(prop))),
                  size=4,colour="white",family = "serif")+
        geom_text(aes(y=rev(rev(prop/2)+c(0,cumsum(rev(prop))[-length(prop)])),
                      x= 1.35,label=process),
                  size=4,color="black",family = "serif")+
        labs(x = '', y = '',title = xlabname[gorder])+
        theme_void()+
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              plot.title  = element_text(hjust = 0.5,family = "serif"))


      pp[[gp]]<-p
    }
    bntp<-patchwork::wrap_plots(pp,nrow = pie.layout.nrow)
  }

  if (geom=="chord") {
    ###chord
    df.plot<-dat[,c(1:2,4)]
    order.sector<- c(unique(df.plot$process)[length(unique(df.plot$process)):1],
                     unique(df.plot$group))

    color_phylum<-color_model
    color_genus<-color_group
    names(color_phylum)<-unique(df.plot$process)
    names(color_genus)<-unique(df.plot$group)
    color.comb<-c(color_phylum,color_genus)
    grid.color<-color.comb[match(order.sector,names(color.comb))]

    circlize::circos.clear()
    par(family="serif",mar=c(2,2,2,2))
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
                            facing = "bending.inside",family = par("serif"),
                            niceFacing = FALSE, adj = c(0.5, -0.5), cex = 1)
    }, bg.border = NA)
    suppressWarnings(warning("circos.track"))
    bntp<-recordPlot()
  }

  if (!is.null(export_path)) {
    filex<-paste(export_path,"/Normalized stochasticity ratio (NST)/NST proportion_",geom,".pdf",sep = "")
    if (geom=="chord") {
      message("Please save the plot mannully!!!")
    } else  ggsave(filex, bntp,height = 9,width = 9)
    cat("\n","Normalized stochasticity ratio (NST) proportion","  has been exported to ","/",filex,"",sep = "","\n")
  } else {
    cat("\n","Normalized stochasticity ratio (NST) proportion","  was not exported to local path","",sep = "","\n")
  }
  return(bntp)
}

