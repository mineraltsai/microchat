"calcMicrochatLefse" <- function(submchat,
                               lda_score=2,
                               file.save=TRUE,
                               export_path="microbial differential analysis") {
  export_path<-paste(export_path,"/microbial differential analysis/LEfSe",sep = "")
  if (file.save) dir.create(export_path, recursive = TRUE)
  if (class(submchat)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  otu_table<-submchat$otu_table
  sel.group<-unique(substr(colnames(otu_table),start = 1,stop = 2))
  tax<-submchat$taxon_table
  sampleda<-data.frame(sample=colnames(otu_table),
                       group=substr(colnames(otu_table),start = 1,stop = 2))
  sampleda<-sampleda%>%tibble::column_to_rownames(var = "sample")

  ps<-phyloseq::phyloseq(
               phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = T),
               phyloseq::tax_table(as.matrix(tax)),
               phyloseq::sample_data(sampleda))

  set.seed(8084)
  deres <- MicrobiotaProcess::diff_analysis(obj=ps, class="group",
                         mlfun="lda",ratio=1,
                         filtermod="pvalue", ##  pvalue
                         firstcomfun = "kruskal.test",  #  oneway.test
                         firstalpha=0.05,
                         strictmod=TRUE,
                         secondcomfun = "wilcox.test",
                         subclmin=3,
                         subclwilc=TRUE,
                         secondalpha=0.05,
                         lda=lda_score)
  deres

  file2=paste(export_path,"/lefse result (lda-",lda_score,") .txt",sep = "")
  if (file.save) write.table(deres@result,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","lefse analysis results (lda: ",lda_score,") has been exported to ","/",export_path,"",sep = "","\n")


  microchatLefseobj<-list(deres,
                             deres@result,
                             deres@kwres,
                             deres@mlres,
                             deres@someparams,
                          all.group,
                          sel.group)

  names(microchatLefseobj)[1]<-"data"
  names(microchatLefseobj)[2]<-"res"
  names(microchatLefseobj)[3]<-"res_diff_first"
  names(microchatLefseobj)[4]<-"res_diff_second"
  names(microchatLefseobj)[5]<-"param"
  names(microchatLefseobj)[6]<-"all.group"
  names(microchatLefseobj)[7]<-"sel.group"

  class(microchatLefseobj) <- c("microchat","data.frame")

  return(microchatLefseobj)
}


"plotMicrochatLefse" <- function(microchatLefseobj,
                                 layout="radial",
                                 color_group=colorCustom(5,pal = "gygn"),
                                 label_plot_display=5,
                                 file.save=TRUE,
                                 export_path="microbial differential analysis") {
  suppressMessages(library(MicrobiotaProcess))
  export_path<-paste(export_path,"/microbial differential analysis/LEfSe",sep = "")
  if (file.save) dir.create(export_path, recursive = TRUE)
  if (class(microchatLefseobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  diffres<-microchatLefseobj$data

  all.group<-microchatLefseobj$all.group
  sel.group<-microchatLefseobj$sel.group

  sel.order<-match(sel.group,all.group)

  color_group<-color_group[sel.order]
  names(color_group)<-sel.group


  p_legend<-MicrobiotaProcess::ggdiffclade(obj=diffres,
                        layout=layout,
                        alpha=0.8,
                        linewd=0.2,
                        skpointsize=0.8,
                        taxlevel=label_plot_display,
                        cladetext=2,
                        setColors=T,
                        removeUnknown=T,
                        reduce=T)+
    scale_fill_manual(values=color_group,name = "\nGroup",
                      guide=guide_legend(order=1,nrow = 1)) +
    scale_size_continuous(range = c(0.5, 2), name = "\nNode size (pvalue)",
                          guide = guide_legend(keywidth = 0.2,
                                               keyheight = 1,
                                               order = 2,nrow = 1))+
    theme(panel.background=element_rect(fill="NA"),
          legend.position="right",
          plot.margin=margin(0,0,0,0),
          axis.title.x=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = -1),
          axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          legend.spacing.y=unit(0.02, "cm"),
          legend.title=element_text(size=12,face = "bold",family = "serif"),
          legend.text=element_text(size=8,face = "bold",family = "serif"),
          legend.box.spacing=unit(0.02,"cm"))

  p_nolegend<-MicrobiotaProcess::ggdiffclade(obj=diffres,
                          layout=layout,
                          alpha=0.8,
                          linewd=0.2,
                          skpointsize=0.8,
                          taxlevel=label_plot_display,
                          cladetext=2,
                          setColors=T,
                          removeUnknown=T,
                          reduce=T)+
    scale_fill_manual(values=color_group,
                      guide=guide_legend(order=1,nrow = 1)) +
    scale_size_continuous(range = c(0.5, 2),
                          guide = guide_legend(keywidth = 0.2,
                                               keyheight = 1,
                                               order = 2,nrow = 1))+
    theme(panel.background=element_rect(fill=NA),
          plot.margin=margin(0,0,0,0),
          axis.title.x=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = -1),
          axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          legend.position="none")

  if (file.save) ggsave(paste(export_path,"/lefse cladogram (with legend)",".pdf",sep = ""),p_legend)
  if (file.save) ggsave(paste(export_path,"/lefse cladogram (without legend)",".pdf",sep = ""),p_nolegend)

  cat("\n","Lefse cladogram of differential analysis has been exported. Please check it.","\n")

  return(list(p_nolegend=p_nolegend,p_legend=p_legend))
}


"plotMicrochatComplexLefse" <- function(submchat,lda_score=2,
                                        control="ct",
                                        layout="radial",
                                        label_plot_display=4,
                                        color_group=c("#FF0000","#333399","#009933","#00FFFF","#EC00FF","yellow"),
                                        export_path="microbial differential analysis") {

  suppressMessages(library(MicrobiotaProcess))
  export_path<-paste(export_path,"/microbial differential analysis/LEfSe/Comparison_",control,sep = "")
  dir.create(export_path, recursive = TRUE)
  if (class(submchat)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  abun<-submchat$otu_table
  sel.group<-unique(substr(colnames(abun),start = 1,stop = 2))
  taxon<-submchat$taxon_table
  tree<-submchat$tree
  param<-submchat$param_table

  groupnames<-substr(colnames(abun),start = 1,stop = 2)%>%unique()

  abun_control<-abun[,which(substr(colnames(abun),start = 1,stop = 2)== control)]

  sel.order<-match(sel.group,all.group)
  color_group<-color_group[sel.order]
  if (length(color_group)<length(groupnames)) message("colors used for coloring group are not enough !!!")
  names(color_group)<-sel.group

  for (kk in setdiff(groupnames,control)) {
    abun_treat<-abun[,which(substr(colnames(abun),start = 1,stop = 2)== kk)]
    abun_new<-cbind(abun_control,abun_treat)
    color_groups<-color_group[match(unique(substr(colnames(abun_new),start = 1,stop = 2)),
                                    names(color_group))]
    mchat.sel<-setMicrochat(abun_new,taxon,tree)

    res<-calcMicrochatLefse(mchat.sel,
                            file.save=FALSE,
                            lda_score=lda_score,
                            export_path=export_path)

    res$res

    file2=paste(export_path, "/lefse (",control,"-",kk,")",".txt",sep = "" )
    write.table(res$res,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
    cat("lefse analysis results have been exported to ","/",export_path,"",sep = "")

    p<-plotMicrochatLefse(res,
                          layout=layout,
                          color_group=color_groups,
                          label_plot_display=label_plot_display,
                          file.save=FALSE,
                          export_path=export_path)

    ggsave(paste(export_path,"/lefse_nolegend (",control,"-",kk,")",".pdf",sep = ""),p$p_nolegend)
    ggsave(paste(export_path,"/lefse_legend (",control,"-",kk,")",".pdf",sep = ""),p$p_legend)
    cat("Pairwise lefse cladogram have been exported to","/",export_path,"",sep = "")
  }
}

"fanfun" <- function(x,pvalue_thres,foldchange_thres){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(pvalue_thres))
  dff <- rbind(data.frame(x = inputx + foldchange_thres, y = y),
               data.frame(x = -(inputx + foldchange_thres), y = y))
  return(dff)
}

"calcMicrochatLimma" <- function(submchat,
                                 comparison = "ml-ct",
                                 file.save=TRUE,
                                 export_path="microbial differential analysis/limma") {
  message(" message: ajust.method=BH,test=glm")
  select_comparison<-strsplit(comparison,"-")
  select_comparison<-select_comparison[[1]]
  message("Control selected is ",select_comparison[1])
  message("Treatment selected is ",select_comparison[2])

  select_comparisonxx<-paste(select_comparison[1],
                           select_comparison[2],sep = "-")

  select_comparison<-paste(select_comparison[2],
                           select_comparison[1],sep = "-")


  export_path<-paste(export_path,"/microbial differential analysis/Limma/",select_comparisonxx,sep = "")
  dir.create(export_path, recursive = TRUE)
  if (class(submchat)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  otu_table<-submchat$otu_table
  tax<-submchat$taxon_table
  sampleda<-data.frame(sample=colnames(otu_table),
                       group=substr(colnames(otu_table),start = 1,stop = 2))
  sampleda<-sampleda%>%tibble::column_to_rownames(var = "sample")

  ps<-phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = T),
    phyloseq::tax_table(as.matrix(tax)),
    phyloseq::sample_data(sampleda))
  counts<-ps@otu_table%>%data.frame()
  group_order<-unique(substr(colnames(counts),start = 1,stop = 2))
  sampledata<-ps@sam_data%>%data.frame()
  sampleda<-rownames_to_column(sampledata,var = "rowname")
  colnames(sampleda)[1]<-"sample"



  my.counts.round<- counts %>% round()

  group_list<-sampleda$group

  edger.counts <- edgeR::DGEList(counts = my.counts.round,
                                 group=factor(group_list),
                                 samples = sampleda)

  filtered.counts <- edger.counts

  filtered.counts <- edgeR::calcNormFactors(filtered.counts, method = c("TMM"))

  design <- model.matrix(~0+factor(filtered.counts$samples$group))
  colnames(design)=levels(factor(filtered.counts$samples$group))
  rownames(design)=colnames(filtered.counts)

  logCPM <- edgeR::cpm(filtered.counts, log=TRUE, prior.count=3)

  v <- limma::voom(filtered.counts,design,plot=FALSE, normalize="quantile")
  fit <- limma::lmFit(v, design)



  cont.matrix=limma::makeContrasts(
    contrasts=select_comparison,
    levels = colnames(design))
  fit2=limma::contrasts.fit(fit,cont.matrix)
  fit2=limma::eBayes(fit2)

  tempOutput = limma::topTable(fit2, coef=select_comparison, n=Inf)
  result = na.omit(tempOutput)
  colnames(result)[2]<-"logCPM"
  colnames(result)[4]<-"pvalue"
  colnames(result)[5]<-"adj_pvalue"
  result<-result%>%tibble::rownames_to_column(var = "name")
  file2=paste(export_path,"/Differential analysis results (limma-- ",select_comparisonxx,") .txt",sep = "")
  if (file.save) write.table(result,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Differential analysis results based on edgeR and limma (",select_comparisonxx,") has been exported to ","/",export_path,"",sep = "","\n")

  microchatLimmaobj<-list(result,
                          select_comparison,
                          submchat$otu_table,
                          submchat$taxon_table)

  names(microchatLimmaobj)[1]<-"result"
  names(microchatLimmaobj)[2]<-"select_comparison"
  names(microchatLimmaobj)[3]<-"otu_table"
  names(microchatLimmaobj)[4]<-"taxon_table"
  class(microchatLimmaobj) <- c("microchat","data.frame")

  return(microchatLimmaobj)
}


"microchatLimmaDiffdata" <- function(microchatLimmaobj,
                                     use.curve=FALSE,
                                     pvalue_thres=0.05,
                                     foldchange_thres=1) {
  data<-microchatLimmaobj
  select_comparison<-data$select_comparison

  control_treat<-strsplit(select_comparison,"-")
  control<-control_treat[[1]][2]
  treat<-control_treat[[1]][1]
  data<-data$result
  data %>%
    mutate(change = case_when(
      pvalue < pvalue_thres & logFC > foldchange_thres ~ "UP",
      pvalue < pvalue_thres & logFC < -foldchange_thres ~ "DOWN",
      TRUE ~ "NOT"
    )) -> DEG

  if (use.curve) {
    if (abs(min(DEG$logFC))>max(DEG$logFC))  xcoordx<-abs(min(DEG$logFC)) else xcoordx<-max(DEG$logFC)
    dff_curve <- fanfun(xcoordx,pvalue_thres,foldchange_thres)
    head(dff_curve)

    DEGx<-DEG
    DEGx$curve_y <- case_when(
      DEGx$logFC > 0 ~ 1/(DEGx$logFC-foldchange_thres) + (-log10(pvalue_thres)),
      DEGx$logFC <= 0 ~ 1/(-DEGx$logFC-foldchange_thres) + (-log10(pvalue_thres))
    )
    DEGx$`-log10(pvalue)`<-(-log10(DEGx$pvalue))
    #根据曲线新增上下调分组标签：
    DEGx$change <- case_when(
      DEGx$`-log10(pvalue)` > DEGx$curve_y & DEGx$logFC >= foldchange_thres ~ 'UP',
      DEGx$`-log10(pvalue)` > DEGx$curve_y & DEGx$logFC <= -foldchange_thres ~ 'DOWN',
      TRUE ~ 'NOT'
    )
    DEG<-DEGx
  }
  DEG$`-log10(pvalue)`<-(-log10(DEG$pvalue))
  DEG$control<-control
  DEG$treat<-treat
  DEG<-DEG[which(DEG$change!="NOT"),]
  return(DEG)
}

"calcDiffData" <- function(submchat,
                           control="CS",
                           pvalue_thres=0.05,
                           foldchange_thres=1,
                           export_path) {
  ggname<-substr(colnames(submchat$otu_table),start = 1,stop = 2)%>%unique()
  diff.data1<-data.frame()
  diff.data2<-data.frame()
  for (i in setdiff(ggname,control)) {
    microchatLimmaobj<-calcMicrochatLimma(submchat,
                                          comparison = paste(control,"-",i,sep = ""),
                                          file.save=TRUE,
                                          export_path=export_path)

    dif.dat1<-microchatLimmaDiffdata(microchatLimmaobj,
                                     use.curve=TRUE,
                                     pvalue_thres=pvalue_thres,
                                     foldchange_thres=foldchange_thres)
    dif.dat2<-microchatLimmaDiffdata(microchatLimmaobj,
                                     use.curve=FALSE,
                                     pvalue_thres=pvalue_thres,
                                     foldchange_thres=foldchange_thres)

    {
      diff.data1<-rbind(diff.data1,dif.dat1)
      diff.data2<-rbind(diff.data2,dif.dat2)
    }
  }

  return(list(diff.data.curve=diff.data1,diff.data.common=diff.data2))
}



"plotMicrochatDiffVolcano" <- function(microchatLimmaobj,
                                       xlabname,
                                       remove.grid=TRUE,
                                       add.barplot=TRUE,
                                       add.pieplot=TRUE,
                                       sigcolor=c('#a121f0','#bebebe','#ffad21'),
                                       sigalpha=c(1,0.4,1),
                                       layout="choice1",
                                       pvalue_thres=0.05,
                                       foldchange_thres=1,
                                       point.size=1,
                                       linetype="dashed",
                                       color_pie.text="white",
                                       color_background="grey90",
                                       export_path="difference analysis/limma"

) {
  export_path<-paste(export_path,"/microbial differential analysis/Limma",sep = "")
  dir.create(export_path, recursive = TRUE)

  if (class(microchatLimmaobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data<-microchatLimmaobj
  select_comparison<-data$select_comparison

  control_treat<-strsplit(select_comparison,"-")
  control<-control_treat[[1]][2]
  treat<-control_treat[[1]][1]

  if (!is.null(names(xlabname))) {
    control<-names(xlabname)[match(control,xlabname)]
    treat<-names(xlabname)[match(treat,xlabname)]
    gp.use<-names(xlabname)
  } else {
    gp.use<-c(control,treat)
  }


  data<-data$result
  data %>%
    mutate(change = case_when(
      pvalue < pvalue_thres & logFC > foldchange_thres ~ "UP",
      pvalue < pvalue_thres & logFC < -foldchange_thres ~ "DOWN",
      TRUE ~ "NOT"
    )) -> DEG
  ylimin<-DEG[which(DEG$change!="NOT"),]$pvalue%>%log10()
  ylimin<-(ylimin*(-1))%>%min()
  label<-table(DEG$change)
  label.o<-c("DOWN","NOT","UP")
  label<-label[match(label.o,names(label))]
  sigcolor<-sigcolor[1:length(unique(names(label)))]
  names(sigcolor)<-unique(names(label))
  names(sigalpha)<-unique(names(label))

  otu2<-DEG[which(DEG$change == "UP"),]
  otu3<-DEG[which(DEG$change == "DOWN"),]
  otu2<-otu2[order(otu2$pvalue),]
  otu3<-otu3[order(otu3$pvalue),]

  otu2<-otu2[1:(as.integer(nrow(otu2)/3)+1),]
  otu3<-otu3[1:(as.integer(nrow(otu3)/3)+1),]

  otu2<-rbind(otu2,otu3)
  if (length(which(is.na(otu2$logFC)))!=0) otu2<-otu2[-which(is.na(otu2$logFC)),]



  if (layout=="choice1") {

    ggplot(data=DEG,aes(x=logFC,
                        y=-log10(pvalue),
                        color=change,alpha=change))+
      geom_point(size=point.size)+
      labs(x="log2 fold change")+ ylab("-log10 pvalue")+
      #ggtitle(this_title)+
      theme_bw(base_size = 10)+
      #theme(plot.title = element_text(size=15,hjust=0.5),)+
      scale_color_manual(values=sigcolor,name="Pattern")+
      scale_alpha_manual(values = sigalpha,guide="none") -> p1

    p2<-p1 +
      geom_vline(xintercept = foldchange_thres,lty=linetype)+
      geom_vline(xintercept = -foldchange_thres,lty=linetype) +
      geom_hline(yintercept = ylimin,lty="dashed")+
      labs(title = paste(treat," regulation patterns compared to ",control,sep = ""))+
      ggrepel::geom_text_repel(data=otu2[which(otu2$logFC>0),],
                               aes(x=logFC,
                                   y=-log10(pvalue), label = name),
                               size = 2.5,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
      ggrepel::geom_text_repel(data=otu2[which(otu2$logFC<0),],
                               aes(x=logFC,
                                   y=-log10(pvalue), label = name),
                               size = 2.5,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))

    p2<-p2+
      annotate('text', size = 4,color=sigcolor[1],
               family="serif",label = paste(treat," Down: ",label[1],sep = ""),
               min(DEG$logFC)*0.8, max((-log10(data$pvalue)))*0.95) +
      annotate('text',size = 4,color=sigcolor[3],
               family="serif", label = paste(treat," Up: ",label[3],sep = ""),
               max(DEG$logFC)*0.9, max((-log10(data$pvalue)))*0.95) +
      theme(aspect.ratio = 1,
            panel.background=element_rect(fill=color_background),
            legend.position = c(.9, .1),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(family = "serif",size=10),
            text = element_text(family = "serif"),
            title = element_text(size = 12),
            axis.text.y=element_text(colour='black',size=10,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 10,
                                     angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
            axis.title.x=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = -1),
            axis.title.y=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5))

    if (remove.grid) p2<-p2+theme(panel.grid = element_blank())

    node_phylum<-table(DEG$change)%>%data.frame()
    node_phylumx<-node_phylum
    node_phylumx<-node_phylumx%>%
      group_by(Var1)%>%
      summarise_all(sum)
    node_phylumx<-node_phylumx[order(node_phylumx$Freq,decreasing = TRUE),]
    node_phylumx$Var1<-factor(node_phylumx$Var1,levels =rev(unique(node_phylumx$Var1)) )

    ylimax<-DEG[which(DEG$change=="NOT"),]$pvalue%>%log10()
    ylimax<-(ylimax*(-1))%>%max()
    yxlim<-max((-log10(data$pvalue)))*0.95/3
    if (ylimax>yxlim) ylimax<-yxlim else ylimax<-ylimax

    ###calc abundance in two group
    sample.size=microchatLimmaobj$otu_table%>%colSums()
    if (length(unique(sample.size))==1) sample.size<-sample.size[1] else sample.size<-max(sample.size)
    otutab<-microchatLimmaobj$otu_table
    split_otux<-calc_indivGroup_abun(otutab)
    split_otux.mean<-(sapply(split_otux, rowMeans)/sample.size*100)%>%data.frame()
    split_otux.mean<-split_otux.mean%>%tibble::rownames_to_column(var = "name")

    abun.sel<-split_otux.mean[,c(1,match(rev(control_treat[[1]]),
                                         colnames(split_otux.mean)))]

    diff.sel<-merge(abun.sel,DEG,by="name")%>%subset(select=c(2,3,10))
    diff.sel$change<-factor(diff.sel$change,levels =label.o )
    plot.data<-diff.sel%>%group_by(change)%>%summarise_all(sum)
    plot.data<-reshape2::melt(plot.data)
    plot.data<-plot.data[order(plot.data$change,
                               plot.data$value,decreasing = TRUE),]
    plot.data$value<-round(plot.data$value,2)
    pp.abun<-ggplot(plot.data,aes(x=variable,y=value))+
      geom_col(position = "stack",show.legend = FALSE,
               aes(fill=change,alpha=change),width = 0.75)+
      geom_text(aes(label=scales::percent(value/sum(value)*2)),
                color=color_pie.text,family="serif",angle=0,
                position = position_stack(vjust = 0.5))+
      scale_fill_manual(values = sigcolor)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      scale_x_discrete(labels = gp.use)+
      theme_void()+
      theme(text = element_text(family = "serif"),
            axis.text.x = element_text(colour='black',
                                       size=12,
                                       family = "serif",vjust = 1),
            legend.position = "none")


    barpv<-ggplot(node_phylum, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', width = 0.5,aes(alpha=Var1)) +
      coord_polar(theta = "y") +
      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.05,
                    label=ifelse((Freq/sum(node_phylum$Freq))>0.04,
                                 scales::percent(Freq/sum(node_phylum$Freq)),"")),
                size=4,colour=color_pie.text,family = "serif")+
      #geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
      #             x= 1.35,label=Var1),
      #         size=4,color="black",family = "serif")+
      scale_fill_manual(values = sigcolor)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank()) +
      labs(x = '', y = '')+theme_void()+
      theme(legend.position = "none")

    x.r<-DEG$logFC%>%max()
    x.l<-DEG$logFC%>%min()
    if (add.pieplot) if (x.r<abs(x.l)) p2<- p2 + annotation_custom(grob=ggplotGrob(barpv),
                                                  ymin = ylimax, ymax = Inf,
                                                  xmin=-Inf,
                                                  xmax=-foldchange_thres*1.1)
    if (add.pieplot) if (x.r>abs(x.l)) p2<- p2 + annotation_custom(grob=ggplotGrob(barpv),
                                                  ymin = ylimax, ymax = Inf,
                                                  xmin=foldchange_thres*1.1,
                                                  xmax=Inf)

    if (add.barplot) p2<- p2 + annotation_custom(grob=ggplotGrob(pp.abun),
                                ymin = ylimax, ymax = Inf,
                                xmin=-foldchange_thres*1.1,
                                xmax=foldchange_thres*1.1)
  }

  if (layout=="choice2"){
    ####choice2
    if (abs(min(DEG$logFC))>max(DEG$logFC))  xcoordx<-abs(min(DEG$logFC)) else xcoordx<-max(DEG$logFC)
    dff_curve <- fanfun(xcoordx,pvalue_thres,foldchange_thres)
    head(dff_curve)

    DEGx<-DEG
    DEGx$curve_y <- case_when(
      DEGx$logFC > 0 ~ 1/(DEGx$logFC-foldchange_thres) + (-log10(pvalue_thres)),
      DEGx$logFC <= 0 ~ 1/(-DEGx$logFC-foldchange_thres) + (-log10(pvalue_thres))
    )
    DEGx$`-log10(pvalue)`<-(-log10(DEGx$pvalue))
    #根据曲线新增上下调分组标签：
    DEGx$change <- case_when(
      DEGx$`-log10(pvalue)` > DEGx$curve_y & DEGx$logFC >= foldchange_thres ~ 'UP',
      DEGx$`-log10(pvalue)` > DEGx$curve_y & DEGx$logFC <= -foldchange_thres ~ 'DOWN',
      TRUE ~ 'NOT'
    )
    labelx<-table(DEGx$change)
    label.o<-c("DOWN","NOT","UP")
    labelx<-labelx[match(label.o,names(labelx))]

    otu2<-DEGx[which(DEGx$change == "UP"),]
    otu3<-DEGx[which(DEGx$change == "DOWN"),]
    otu2<-otu2[order(otu2$pvalue),]
    otu3<-otu3[order(otu3$pvalue),]

    otu2<-otu2[1:(as.integer(nrow(otu2)/3)+1),]
    otu3<-otu3[1:(as.integer(nrow(otu3)/3)+1),]

    otu2<-rbind(otu2,otu3)

    p1 <- ggplot(data = DEGx,
                 aes(x = logFC,
                     y = -log10(pvalue),
                     color = change)) +
      geom_point(size = point.size,aes(alpha=change)) +
      scale_y_continuous(limits = c(NA, max(DEGx$`-log10(pvalue)`)*1.2)) +
      scale_color_manual(values=sigcolor,name="Pattern")+
      scale_alpha_manual(values = sigalpha,guide="none")+
      labs(x="log2 fold change")+ ylab("-log10 pvalue")+
      theme_bw(base_size = 10)+
      geom_line(data = dff_curve,
                aes(x = x, y = y),
                color = "black",lty = linetype)+
      labs(title = paste(treat," regulation patterns compared to ",control,sep = ""))+
      ggrepel::geom_text_repel(data=otu2[which(otu2$logFC>0),],
                               aes(x=logFC,
                                   y=-log10(pvalue), label = name),
                               size = 2.5,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
      ggrepel::geom_text_repel(data=otu2[which(otu2$logFC<0),],
                               aes(x=logFC,
                                   y=-log10(pvalue), label = name),
                               size = 2.5,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))


    if (abs(min(DEGx$logFC))>max(DEGx$logFC))  xcoord<-abs(min(DEGx$logFC)) else xcoord<-max(DEGx$logFC)
    p2<-p1+
      annotate('text', size = 4,color=sigcolor[1],
               family="serif",label = paste(treat," Down: ",labelx[1],sep = ""),
               -xcoord*0.9, max((-log10(data$pvalue)))*1.1) +
      annotate('text',size = 4,color=sigcolor[3],
               family="serif", label = paste(treat," Up: ",labelx[3],sep = ""),
               xcoord*0.9, max((-log10(data$pvalue)))*1.1) +
      theme(aspect.ratio = 1,
            panel.background=element_rect(fill=color_background),
            legend.position = c(.9, .1),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(family = "serif",size=10),
            text = element_text(family = "serif"),
            title = element_text(size = 12),
            axis.text.y=element_text(colour='black',size=10,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 10,
                                     angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
            axis.title.x=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = -1),
            axis.title.y=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5))


    if (remove.grid) p2<-p2+theme(panel.grid = element_blank())

    node_phylum<-table(DEGx$change)%>%data.frame()
    node_phylumx<-node_phylum
    node_phylumx<-node_phylumx%>%
      group_by(Var1)%>%
      summarise_all(sum)
    node_phylumx<-node_phylumx[order(node_phylumx$Freq,decreasing = TRUE),]
    node_phylumx$Var1<-factor(node_phylumx$Var1,levels =rev(unique(node_phylumx$Var1)) )

    ylimax<-DEGx[which(DEGx$change=="NOT"),]$pvalue%>%log10()
    ylimax<-(ylimax*(-1))%>%max()
    yxlim<-max((-log10(data$pvalue)))*0.95/3
    if (ylimax>yxlim) ylimax<-yxlim else ylimax<-ylimax

    ###calc abundance in two group
    sample.size=microchatLimmaobj$otu_table%>%colSums()
    if (length(unique(sample.size))==1) sample.size<-sample.size[1] else sample.size<-max(sample.size)
    otutab<-microchatLimmaobj$otu_table
    split_otux<-calc_indivGroup_abun(otutab)
    split_otux.mean<-(sapply(split_otux, rowMeans)/sample.size*100)%>%data.frame()
    split_otux.mean<-split_otux.mean%>%tibble::rownames_to_column(var = "name")

    abun.sel<-split_otux.mean[,c(1,match(rev(control_treat[[1]]),
                                         colnames(split_otux.mean)))]

    diff.sel<-merge(abun.sel,DEGx,by="name")%>%subset(select=c(2,3,10))
    diff.sel$change<-factor(diff.sel$change,levels =label.o )
    plot.data<-diff.sel%>%group_by(change)%>%summarise_all(sum)
    plot.data<-reshape2::melt(plot.data)
    plot.data<-plot.data[order(plot.data$change,
                               plot.data$value,decreasing = TRUE),]
    plot.data$value<-round(plot.data$value,2)
    pp.abun<-ggplot(plot.data,aes(x=variable,y=value))+
      geom_col(position = "stack",show.legend = FALSE,
               aes(fill=change,alpha=change),width = 0.75)+
      geom_text(aes(label=scales::percent(value/sum(value)*2)),
                color=color_pie.text,family="serif",angle=0,
                position = position_stack(vjust = 0.5))+
      scale_fill_manual(values = sigcolor)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      scale_x_discrete(labels = gp.use)+
      theme_void()+
      theme(text = element_text(family = "serif"),
            axis.text.x = element_text(colour='black',
                                       size=12,
                                       family = "serif",vjust = 1),
            legend.position = "none")


    barpv<-ggplot(node_phylum, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', width = 0.5,aes(alpha=Var1)) +
      coord_polar(theta = "y") +
      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.05,
                    label=ifelse((Freq/sum(node_phylum$Freq))>0.04,
                                 scales::percent(Freq/sum(node_phylum$Freq)),"")),
                size=4,colour=color_pie.text,family = "serif")+
      # geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
      #              x= 1.35,label=Var1),
      #          size=4,color="black",family = "serif")+
      scale_fill_manual(values = sigcolor)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_blank()) +
      labs(x = '', y = '')+theme_void()+
      theme(legend.position = "none")

    if (add.pieplot)  p2<- p2 + annotation_custom(grob=ggplotGrob(barpv),
                                                  ymin = ylimax, ymax = Inf,
                                                  xmin=foldchange_thres*1.1,
                                                  xmax=Inf)

  if (add.barplot)  p2<- p2 + annotation_custom(grob=ggplotGrob(pp.abun),
                                ymin = ylimax,
                                ymax = max((-log10(data$pvalue)))*1.1,
                                xmin=-Inf,
                                xmax=-foldchange_thres*2)
   DEG<-DEGx
  }

  print(p2)

  ggsave(paste(export_path,"/single_volcano_",layout,"_",paste(treat,control,sep = "-"),".pdf",sep = ""),
         width = 6,height = 6,
         p2)
  cat("Single volcano has been exported to","/",export_path,"",sep = "")

  taxo<-microchatLimmaobj$taxon_table%>%tibble::rownames_to_column(var = "name")
  DEG<-merge(DEG,abun.sel,by="name")
  DEG<-merge(DEG,taxo,by="name")
  #DEG$control<-control
  #DEG$treat<-treat
  file2=paste(export_path, "/single_volcano_",paste(control,treat,sep = "-"),"_oridata.txt",sep = "" )
  write.table(DEG,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  file2=paste(export_path, "/single_volcano_",paste(control,treat,sep = "-"),"_onlydiffdata.txt",sep = "" )
  write.table(DEG[which(DEG$change!="NOT"),],file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  cat("\n","Differential analysis results has been exported to“",export_path,"",sep = "")

  return(list(data=DEG,plot=p2,sig_condtion=label,diff.data=DEG[which(DEG$change!="NOT"),]))
}


"edger_mutigroup" <- function(counts,sampleda,
                              select_comparison) {

  library(tidyverse)
  library(edgeR)
  my.counts.round<- counts %>% round()

  group_list<-sampleda$group

  edger.counts <- DGEList(counts = my.counts.round,
                          group=factor(group_list),
                          samples = sampleda)

  filtered.counts <- edger.counts

  filtered.counts <- calcNormFactors(filtered.counts, method = c("TMM"))

  design <- model.matrix(~0+factor(filtered.counts$samples$group))
  colnames(design)=levels(factor(filtered.counts$samples$group))
  rownames(design)=colnames(filtered.counts)

  logCPM <- cpm(filtered.counts, log=TRUE, prior.count=3)

  v <- voom(filtered.counts,design,plot=FALSE, normalize="quantile")
  fit <- lmFit(v, design)

  select_comparison<-paste(select_comparison[2],
                           select_comparison[1],sep = "-")

  cont.matrix=makeContrasts(
    contrasts=select_comparison,
    levels = colnames(design))
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)

  tempOutput = topTable(fit2, coef=select_comparison, n=Inf)
  result = na.omit(tempOutput)
  colnames(result)[2]<-"logCPM"
  colnames(result)[4]<-"pvalue"
  colnames(result)[5]<-"adj_pvalue"


  return(result)
}



"calcMicrochatComplexLimma" <- function(submchat,
                                        comparison = "ct-xl-yl",
                                        file.save=TRUE,
                                        export_path="microbial differential analysis/limma") {
  message(" message: ajust.method=BH,test=glm")
  export_path<-paste(export_path,"/microbial differential analysis/Limma",sep = "")
  dir.create(export_path, recursive = TRUE)
  if (class(submchat)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  otu_table<-submchat$otu_table
  tax<-submchat$taxon_table
  sampleda<-data.frame(sample=colnames(otu_table),
                       group=substr(colnames(otu_table),start = 1,stop = 2))
  sampleda<-sampleda%>%tibble::column_to_rownames(var = "sample")

  ps<-phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = T),
    phyloseq::tax_table(as.matrix(tax)),
    phyloseq::sample_data(sampleda))
  counts<-ps@otu_table%>%data.frame()
  group_order<-unique(substr(colnames(counts),start = 1,stop = 2))
  sampledata<-ps@sam_data%>%data.frame()
  sampleda<-rownames_to_column(sampledata,var = "rowname")
  colnames(sampleda)[1]<-"sample"


  ###chose comparison
  select_comparison<-strsplit(comparison,"-")

  select_comparison<-select_comparison[[1]]

  select_comparison1<-select_comparison[1:2]
  select_comparison2<-cbind(select_comparison[1],select_comparison[3])%>%as.character()


  comparison11<-edger_mutigroup(counts=counts,
                                sampleda=sampleda,
                                select_comparison=select_comparison1)

  comparison21<-edger_mutigroup(counts=counts,
                                sampleda=sampleda,
                                select_comparison=select_comparison2)

  comparison11<-subset(comparison11,select = -c(t,B))
  comparison21<-subset(comparison21,select = -c(t,B))

  colnames(comparison11)<-paste(colnames(comparison11),"_1",sep = "")
  colnames(comparison21)<-paste(colnames(comparison21),"_2",sep = "")

  comparison11<-rownames_to_column(comparison11,var = "rownames")
  comparison21<-rownames_to_column(comparison21,var = "rownames")

  cpma<-merge(comparison11,comparison21,by="rownames")
  colnames(cpma)
  cpma<-subset(cpma,select = -logCPM_2)
  colnames(cpma)[3]<-"logCPM"


  file2=paste(export_path,"/Muti-differential analysis results (limma-- ",comparison,") .txt",sep = "")
  if (file.save) write.table(cpma,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Muti-differential analysis results based on edgeR and limma (",select_comparison,") has been exported to ","/",export_path,"",sep = "","\n")

  microchatComplexLimmaobj<-list(data=cpma,
                                 select_comparison1=select_comparison1,
                                 select_comparison2=select_comparison2,
                                 sampleda=sampleda,
                                 otu_table=otu_table)


  names(microchatComplexLimmaobj)[1]<-"data"
  names(microchatComplexLimmaobj)[2]<-"select_comparison1"
  names(microchatComplexLimmaobj)[3]<-"select_comparison2"
  names(microchatComplexLimmaobj)[4]<-"sampleda"
  names(microchatComplexLimmaobj)[5]<-"otu_table"

  class(microchatComplexLimmaobj) <- c("microchat","data.frame")

  return(microchatComplexLimmaobj)
}


"plotMicrochatComplexDiffVolcano" <- function(microchatComplexLimmaobj,
                                              add.barplot=TRUE,
                                              repelsize=3,
                                              color_pie.text="white",
                                              repeldisp_rate=2, #1-all,2-half
                                              annosize=3,
                                              pvalue=0.05,
                                              foldchange=1,
                                              sigcolor=c('red', 'orange', 'green', 'blue'),
                                              export_path="difference analysis/limma") {
  export_path<-paste(export_path,"/microbial differential analysis/Limma",sep = "")
  dir.create(export_path, recursive = TRUE)

  if (class(microchatComplexLimmaobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data<-microchatComplexLimmaobj

  select_comparison1=data$select_comparison1
  select_comparison2=data$select_comparison2
  data<-data$data
  control<-select_comparison1[1]
  treat1<-select_comparison1[2]
  treat2<-select_comparison2[2]
  data[which(data$adj_pvalue_1 < pvalue & data$adj_pvalue_2 < pvalue),'type1'] <- 'sign'
  data[which(data$adj_pvalue_1 >= pvalue | data$adj_pvalue_2 >= pvalue),'type1'] <- 'no'

  #1
  data[which(data$logFC_1 <= -foldchange & data$logFC_2 <= -foldchange),'type2'] <- '1_down.2_down'
  data[which(data$logFC_1 >= foldchange & data$logFC_2 <= -foldchange),'type2'] <- '1_up.2_down'
  data[which(data$logFC_1 <= -foldchange & data$logFC_2 >= foldchange),'type2'] <- '1_down.2_up'
  data[which(data$logFC_1 >= foldchange & data$logFC_2 >= foldchange),'type2'] <- '1_up_2_up'
  data[is.na(data$type2),'type2'] <- 'no'

  #2
  data$type3 <- paste(data$type1, data$type2, sep = '.')

  #3
  data$type3 <- factor(data$type3, levels = c('sign.1_down.2_down', 'sign.1_up.2_down', 'sign.1_down.2_up', 'sign.1_up_2_up',
                                              'no.1_down.2_down', 'no.1_up.2_down', 'no.1_down.2_up', 'no.1_up_2_up',
                                              'sign.no', 'no.no'))
  data <- data[order(data$type3, decreasing = TRUE), ]
  #4
  p <- ggplot(data, aes(logFC_1, logFC_2)) +
    geom_point(aes(color = type3, size = logCPM), alpha = 0.6) +
    scale_size(range = c(0, 4)) +
    scale_color_manual(limits = c('sign.1_down.2_down', 'sign.1_up.2_down', 'sign.1_down.2_up', 'sign.1_up_2_up',
                                  'no.1_down.2_down', 'no.1_up.2_down', 'no.1_down.2_up', 'no.1_up_2_up',
                                  'sign.no', 'no.no'),
                       values = c(sigcolor, 'gray', 'gray', 'gray', 'gray', 'gray', 'gray'),
                       guide=guide_legend(order=1,ncol = 2)) +
    guides(color="none")+
    theme(aspect.ratio = 1,
          text = element_text(family = "serif"),
          axis.text.y=element_text(colour='black',size=10,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 10,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          axis.title.x=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = -1),
          axis.title.y=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
          panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    geom_vline(xintercept = c(-foldchange, foldchange), lty = 2) +
    geom_hline(yintercept = c(-foldchange, foldchange), lty = 2) +
    labs(x = paste('log2FC (',treat1 ,'vs',control,')'), y = paste('log2FC (',treat2 ,'vs',control,')'))


  #5
  p1 <- p +
    annotate('text', size = annosize,family="serif",label = paste(treat1,' Down\n',treat2,' Down',sep = ""), min(data$logFC_1)*1.1, min(data$logFC_2)*1.1) +
    annotate('text',size = annosize,family="serif", label = paste(treat1,' Down\n',treat2,' Up',sep = ""), min(data$logFC_1)*1.1, max(data$logFC_2)*1.1) +
    annotate('text',size = annosize, family="serif",label = paste(treat1,' Up\n',treat2,' Down',sep = ""), max(data$logFC_1)*1.1, min(data$logFC_2)*1.1) +
    annotate('text', size = annosize,family="serif",label = paste(treat1,' Up\n',treat2,' Up',sep = ""), max(data$logFC_1)*1.1, max(data$logFC_2)*1.1)

  data[which(data$adj_pvalue_1 < pvalue & data$logFC_1 <= -foldchange),'sig'] <- 'Down'
  data[which(data$adj_pvalue_1 < pvalue & data$logFC_1 >= foldchange),'sig'] <- 'Up'
  data[which(data$adj_pvalue_2 < pvalue & data$logFC_2 <= -foldchange),'sig'] <- 'Down'
  data[which(data$adj_pvalue_2 < pvalue & data$logFC_2 >= foldchange),'sig'] <- 'Up'
  data[which(data$adj_pvalue_1 >= pvalue | abs(data$logFC_1) < foldchange),'sig'] <- 'None'
  data[which(data$adj_pvalue_2 >= pvalue | abs(data$logFC_2) < foldchange),'sig'] <- 'None'

  up <- subset(data, sig == 'Up')
  up_num<-(length(up$rownames)/repeldisp_rate)%>%as.integer()
  up <- up[order(up$adj_pvalue_1,up$adj_pvalue_2), ][1:up_num, ]

  down <- subset(data, sig == 'Down')
  down_num<-(length(down$rownames)/repeldisp_rate)%>%as.integer()
  down <- down[order(down$adj_pvalue_1,down$adj_pvalue_2), ][1:down_num, ]

  #6
 lgd.x<-1-(max(data$logFC_1)+foldchange)/(max(data$logFC_1)-min(data$logFC_1))
 lgd.x<-round(lgd.x*0.8,2)
 p2 <- p1 +
    theme(legend.position = c(lgd.x,0.1),
          legend.background = element_blank(),
          legend.key = element_blank()) +
    ggrepel::geom_text_repel(data = rbind(up, down),
                             aes(x = logFC_1, y = logFC_2, label = rownames),
                             size = repelsize,box.padding = unit(0.5, 'lines'), family="serif",
                             segment.color = 'black', show.legend = FALSE)

 diff.data<-data[which(data$type1=="sign"),]
 colnames(diff.data)[1]<-"name"

 sample.size=microchatComplexLimmaobj$otu_table%>%colSums()
 if (length(unique(sample.size))==1) sample.size<-sample.size[1] else sample.size<-max(sample.size)
 otutab<-microchatComplexLimmaobj$otu_table
 split_otux<-calc_indivGroup_abun(otutab)
 split_otux.mean<-(sapply(split_otux, rowMeans))%>%data.frame()
 split_otux.mean<-split_otux.mean%>%tibble::rownames_to_column(var = "name")

 abun.sel<-split_otux.mean[,c(1,match(c(control,treat1,treat2),
                                      colnames(split_otux.mean)))]

 diff.sel<-merge(abun.sel,diff.data,by="name")%>%subset(select=c(2,3,4,14))

 sigcolorx<-sigcolor
 names(sigcolorx)<-c('sign.1_down.2_down',
                    'sign.1_up.2_down',
                    'sign.1_down.2_up',
                    'sign.1_up_2_up')
 diff.sel$type3<-factor(diff.sel$type3,levels =unique(diff.sel$type3) )
 plot.data<-diff.sel%>%group_by(type3)%>%summarise_all(sum)
 plot.data<-reshape2::melt(plot.data)
 plot.data<-plot.data[order(plot.data$type3,
                            plot.data$value,decreasing = TRUE),]
 plot.data$value<-round(plot.data$value,2)
 plot.data<-plot.data[-which(plot.data$type3=="sign.no"),]
 plot.data<-plot.data%>%group_by(variable)%>%mutate(prop=value/sum(value)*100)

 plot.data$prop<-round(plot.data$prop,2)
 pp.abun<-ggplot(plot.data,aes(x=variable,y=prop))+
   geom_col(position = "stack",show.legend = FALSE,
            aes(fill=type3),width = 0.75)+
   geom_text(aes(label=paste(prop,"%",sep="")),
             color=color_pie.text,family="serif",angle=0,
             position = position_stack(vjust = 0.5))+
   scale_fill_manual(values = sigcolorx)+
   theme_void()+
   theme(text = element_text(family = "serif"),
         axis.text.x = element_text(colour='black',
                                    size=12,
                                    family = "serif",vjust = 1),
         legend.position = "none")

 if (add.barplot) p2<-p2+annotation_custom(grob=ggplotGrob(pp.abun),
                          ymin = -Inf, ymax =-foldchange,
                          xmin=foldchange,
                          xmax=max(data$logFC_1))


  ggsave(paste(export_path,"/Muti-volcano (",control,"-",treat1,"-",treat2,").pdf",sep = ""),p2)
  cat("Muti-volcano plot has been exported to","/",export_path,"",sep = "")
  return(p2)

}


"plotMicrochatDiffCombinedVolcano" <- function(submchat,
                                             ncol=4,
                                             my_comparisons,
                                             xlabname=NULL,
                                             remove.grid=TRUE,
                                             add.barplot=TRUE,
                                             sigcolor=colorCustom(3,pal = "set1"),
                                             sigalpha=c(1,0.3,1),
                                             layout="choice2",
                                             pvalue_thres=0.05,
                                             foldchange_thres=1,
                                             point.size=1,
                                             linetype="dashed",
                                             color_pie.text="black",
                                             color_background=NA,
                                             export_path) {

  sample.dat<-group_generate(submchat$otu_table)
  g.name<-sample.dat$group%>%unique()
  names(g.name)<-xlabname

  xp<-list()
  for (i in 1:length(my_comparisons)) {

    g.name.sel<-my_comparisons[[i]]
    xlabname.sel<-g.name[match(g.name.sel,g.name)]

    microchatLimmaobj<-calcMicrochatLimma(submchat,
                                          comparison = paste(my_comparisons[[i]][1],
                                                             "-",
                                                             my_comparisons[[i]][2],sep = ""),
                                          file.save=TRUE,
                                          export_path=export_path)

    ### single volcano
    ptreat2<-plotMicrochatDiffVolcano(microchatLimmaobj,
                                      xlabname=xlabname.sel,
                                      remove.grid=remove.grid,
                                      add.barplot=add.barplot,
                                      sigcolor=sigcolor,
                                      sigalpha=sigalpha,
                                      layout=layout,
                                      pvalue_thres=pvalue_thres,
                                      foldchange_thres=foldchange_thres,
                                      point.size=point.size,
                                      linetype=linetype,
                                      color_pie.text=color_pie.text,
                                      color_background=color_background,
                                      export_path=export_path)

    xp[[i]]<-ptreat2$plot
  }

  if (!is.null(ncol)){
    width.sel<-ncol
    if ((length(xp)/width.sel)%%1==0) {
      height.sel<-(length(xp)/width.sel)%>%as.integer()
    } else {
      height.sel<-(length(xp)/width.sel)%>%as.integer()+1
    }
  } else {
    width.sel<-ncol_layout(length(xp))[[1]]%>%length()
    height.sel<-ncol_layout(length(xp))%>%length()
  }

  mp<-patchwork::wrap_plots(xp,ncol = width.sel)
  ggsave(paste(export_path,"/microbial differential analysis/Limma/Combined_volcano.tiff",
               sep = ""),
         width = 6*width.sel,height = 6*height.sel,
         mp)
  ggsave(paste(export_path,"/microbial differential analysis/Limma/Combined_volcano.pdf",
               sep = ""),
         width = 6*width.sel,height = 6*height.sel,
         mp)

}
