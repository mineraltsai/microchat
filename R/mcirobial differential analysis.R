"calcMicrochatLefse" <- function(submchat,
                               lda_score=2,
                               file.save=TRUE,
                               export_path="microbial differential analysis") {
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/LEfSe",sep = "")
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
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/LEfSe",sep = "")
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

  #cat("\n","Lefse cladogram of differential analysis has been exported. Please check it.","\n")

  return(list(p_nolegend=p_nolegend,p_legend=p_legend))
}


"plotMicrochatComplexLefse" <- function(submchat,lda_score=2,
                                        control="ct",
                                        layout="radial",
                                        label_plot_display=4,
                                        color_group=c("#FF0000","#333399","#009933","#00FFFF","#EC00FF","yellow"),
                                        export_path="microbial differential analysis") {

  suppressMessages(library(MicrobiotaProcess))
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/LEfSe/Comparison_",control,sep = "")
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
    #cat("lefse analysis results have been exported to ","/",export_path,"",sep = "")

    p<-plotMicrochatLefse(res,
                          layout=layout,
                          color_group=color_groups,
                          label_plot_display=label_plot_display,
                          file.save=FALSE,
                          export_path=export_path)

    ggsave(paste(export_path,"/lefse_nolegend (",control,"-",kk,")",".pdf",sep = ""),p$p_nolegend)
    ggsave(paste(export_path,"/lefse_legend (",control,"-",kk,")",".pdf",sep = ""),p$p_legend)
    #cat("Pairwise lefse cladogram have been exported to","/",export_path,"",sep = "")
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


  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/Limma/",select_comparisonxx,sep = "")
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
  #cat("\n","Differential analysis results based on edgeR and limma (",select_comparisonxx,") has been exported to ","/",export_path,"",sep = "","\n")

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
                                       color_group=NULL,
                                       sigalpha=c(1,0.4,1),
                                       layout="choice1",
                                       pvalue_thres=0.05,
                                       foldchange_thres=1,
                                       point.size=1,
                                       linetype="dashed",
                                       color_pie.text="white",
                                       color_background="grey90",
                                       mytheme=NULL,
                                       export_path="difference analysis/limma"

) {
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/Limma",sep = "")
  dir.create(export_path, recursive = TRUE)

  if (class(microchatLimmaobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data<-microchatLimmaobj
  select_comparison<-data$select_comparison
  control_treat<-strsplit(select_comparison,"-")
  control<-control_treat[[1]][2]
  treat<-control_treat[[1]][1]

  gxname<-group_generate(data$otu_table)$group%>%unique()

  if (!is.null(color_group)) sigcolor<-c(color_group[match(control,gxname)],"grey75",color_group[match(treat,gxname)]) else sigcolor<-sigcolor

  if (!is.null(xlabname)) {
    names(xlabname)<-group_generate(data$otu_table)$group%>%unique()
    control<-xlabname[match(control,names(xlabname))]
    treat<-xlabname[match(treat,names(xlabname))]
    gp.use<-xlabname
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
      geom_point(data = DEG[which(DEG$change == "NOT"),],
                 size = point.size,
                 aes(alpha=change)) +
      geom_point(data = DEG[which(DEG$change != "NOT"),],
                 size = point.size*1.5,
                 aes(alpha=change)) +
      labs(x="log2 fold change")+ ylab("-log10 pvalue")+
      #ggtitle(this_title)+
      theme_bw(base_size = 10)+
      #theme(plot.title = element_text(size=15,hjust=0.5),)+
      scale_color_manual(values=sigcolor,name="Pattern",guide="none")+
      scale_alpha_manual(values = sigalpha,guide="none") -> p1

    p2<-p1 +
      geom_vline(xintercept = foldchange_thres,lty=linetype)+
      geom_vline(xintercept = -foldchange_thres,lty=linetype) +
      geom_hline(yintercept = ylimin,lty="dashed")+
      labs(title = paste(treat," regulation patterns compared to ",control,sep = ""))+
      ggrepel::geom_text_repel(data=otu2[which(otu2$logFC>0),],
                               aes(x=logFC,
                                   y=-log10(pvalue), label = name),
                               size = 2,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
      ggrepel::geom_text_repel(data=otu2[which(otu2$logFC<0),],
                               aes(x=logFC,
                                   y=-log10(pvalue), label = name),
                               size = 2,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))

    p2<-p2+
      annotate('text', size = 2.5,color=sigcolor[1],fontface = "bold",
               family="serif",label = paste(treat," Down: ",label[1],sep = ""),
               min(DEG$logFC)*0.8, max((-log10(data$pvalue)))*0.98) +
      annotate('text',size = 2.5,color=sigcolor[3],fontface = "bold",
               family="serif", label = paste(treat," Up: ",label[3],sep = ""),
               max(DEG$logFC)*0.9, max((-log10(data$pvalue)))*0.98) +
      theme(aspect.ratio = 1,
            panel.background=element_rect(fill=color_background),
            legend.position = "none",
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
      geom_text(aes(label=paste0(sprintf("%0.2f", value/sum(value)*2*100),"%")),size=2,
                color=color_pie.text,family="serif",angle=0,
                position = position_stack(vjust = 0.5))+
      scale_fill_manual(values = sigcolor)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      scale_x_discrete(labels = gp.use)+
      theme_void()+
      theme(text = element_text(family = "serif"),
            axis.text.x = element_text(colour='black',
                                       size=6,face = "bold",
                                       family = "serif",vjust = 1),
            legend.position = "none")


    barpv<-ggplot(node_phylum, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', width = 0.5,aes(alpha=Var1)) +
      coord_polar(theta = "y") +
      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.05,
                    label=ifelse((Freq/sum(node_phylum$Freq))>0.04,
                                 paste0(sprintf("%0.2f", Freq/sum(node_phylum$Freq)*100),"%"),"")),
                size=2,colour=color_pie.text,family = "serif")+
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
      geom_point(data = DEGx[which(DEGx$change == "NOT"),],
                 size = point.size,
                 aes(alpha=change)) +
      geom_point(data = DEGx[which(DEGx$change != "NOT"),]
                 ,size = point.size*1.5,
                 aes(alpha=change)) +
      scale_y_continuous(limits = c(NA, max(DEGx$`-log10(pvalue)`)*1.2)) +
      scale_color_manual(values=sigcolor,name="Pattern",guide="none")+
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
                               size = 2,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
      ggrepel::geom_text_repel(data=otu2[which(otu2$logFC<0),],
                               aes(x=logFC,
                                   y=-log10(pvalue), label = name),
                               size = 2,family="serif",
                               box.padding = unit(0.3, 'lines'),
                               show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))


    if (abs(min(DEGx$logFC))>max(DEGx$logFC))  xcoord<-abs(min(DEGx$logFC)) else xcoord<-max(DEGx$logFC)
    p2<-p1+
      annotate('text', size = 2.5,color=sigcolor[1],fontface = "bold",
               family="serif",label = paste(treat," Down: ",labelx[1],sep = ""),
               -xcoord*0.9, max((-log10(data$pvalue)))*1.2) +
      annotate('text',size = 2.5,color=sigcolor[3],fontface = "bold",
               family="serif", label = paste(treat," Up: ",labelx[3],sep = ""),
               xcoord*0.9, max((-log10(data$pvalue)))*1.2) +
      theme(aspect.ratio = 1,
            panel.background=element_rect(fill=color_background),
            legend.position = "none",
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
      geom_text(aes(label=paste0(sprintf("%0.2f", value/sum(value)*2*100),"%")),size=2,
                color=color_pie.text,family="serif",angle=0,
                position = position_stack(vjust = 0.5))+
      scale_fill_manual(values = sigcolor)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      scale_x_discrete(labels = gp.use)+
      theme_void()+
      theme(text = element_text(family = "serif"),
            axis.text.x = element_text(colour='black',
                                       size=6,face = "bold",
                                       family = "serif",vjust = 1),
            legend.position = "none")


    barpv<-ggplot(node_phylum, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(stat = 'identity', width = 0.5,aes(alpha=Var1)) +
      coord_polar(theta = "y") +
      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.05,
                    label=ifelse((Freq/sum(node_phylum$Freq))>0.04,
                                 paste0(sprintf("%0.2f", Freq/sum(node_phylum$Freq)*100),"%"),"")),
                size=2,colour=color_pie.text,family = "serif")+
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

  p2<-p2+theme(aspect.ratio = 1,
             title =  element_text(size = 7.5,face="bold"),
             axis.title = element_text(size = 6,face="bold"),
             axis.text = element_text(size = 4,face="bold"),
             panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5))
  if (!is.null(mytheme)) p2<-p2+mytheme

  taxo<-microchatLimmaobj$taxon_table%>%tibble::rownames_to_column(var = "name")
  DEG<-merge(DEG,abun.sel,by="name")
  DEG<-merge(DEG,taxo,by="name")

  #p2<-p2+labs(tag = paste0("Selected thres\npvalue:"," ",pvalue_thres,"\nfold change: ",foldchange_thres))
  print(p2)
  ggsave(paste(export_path,"/single_volcano_",layout,"_",paste(treat,control,sep = "-"),"_pvalue(",pvalue_thres,")_fc(",sprintf("%0.2f",foldchange_thres),").pdf",sep = ""),
         units = "cm",
         width = 21/3,
         height = 21*p2$theme$aspect.ratio/3,
         p2)
  #cat("Single volcano has been exported to","/",export_path,"",sep = "")

  #DEG$control<-control
  #DEG$treat<-treat
  file2=paste(export_path, "/single_volcano_",paste(control,treat,sep = "-"),"_pvalue(",pvalue_thres,")_fc(",sprintf("%0.2f",foldchange_thres),")_oridata.txt",sep = "" )
  write.table(DEG,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  file2=paste(export_path, "/single_volcano_",paste(control,treat,sep = "-"),"_pvalue(",pvalue_thres,")_fc(",sprintf("%0.2f",foldchange_thres),")_onlydiffdata.txt",sep = "" )
  write.table(DEG[which(DEG$change!="NOT"),],file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  #cat("\n","Differential analysis results has been exported to“",export_path,"",sep = "")

  return(list(data=DEG,plot=p2,microchatLimmaobj=microchatLimmaobj,sigcolor=sigcolor,treat=treat,control=control,sig_condtion=label,pvalue_thres=pvalue_thres,foldchange_thres=foldchange_thres,diff.data=DEG[which(DEG$change!="NOT"),]))
}


'plotMicrochatDiffTaxaDistri' <- function(ptreat,
                                          taxa='Phylum',
                                          sigalpha=c(0.04,0.04),
                                          color_taxa=colorCustom(50,pal = "gygn"),
                                          dwn_coef=0.4,
                                          up_coef=1.35,
                                          legend.position = c(0.35,0.7),
                                          clip_layout=FALSE,  ##是否旋转下层节点分布图
                                          export_path) {
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/Limma",sep = "")
  dir.create(export_path, recursive = TRUE)

  DEG<-ptreat$data
  microchatLimmaobj<-ptreat$microchatLimmaobj
  control<-ptreat$control
    treat<-ptreat$treat
    sigcolor<-ptreat$sigcolor

  DEG[,taxa][which(DEG$change=="UP")]%>%table()%>%as.data.frame()->d1d
  DEG[,taxa][which(DEG$change=="DOWN")]%>%table()%>%as.data.frame()->d2d

  if (nrow(d2d)!=0 | nrow(d1d)!=0) {
    if (nrow(d1d)==0 & nrow(d2d)!=0) {
      d1d<-d2d
      d1d$Freq<-0
    }
    if (nrow(d2d)==0 & nrow(d1d)!=0) {
      d2d<-d1d
      d2d$Freq<-0
    }

    d2d$Freq<-d2d$Freq*(-1)
    colnames(d2d)<-colnames(d1d)<-c("tax",'value')
    d1d$variable<-"UP";d2d$variable<-"UP"
    split_otu4<-rbind(d1d,d2d)

    if ((nrow(split_otu4)-1)>=15) {
      split_otu4<-rbind(data.frame(tax=c("___AAA","___AAB"),
                                   value=c(0,0),
                                   variable=c("UP","UP")),
                        split_otu4)
      split_otu4$label<-abs(split_otu4$value)
      color_taxa_use<-c(NA,NA,color_taxa[1:length(unique(microchatLimmaobj$taxon_table[,taxa]))])
      names(color_taxa_use)<-c("___AAA","___AAB",unique(microchatLimmaobj$taxon_table[,taxa]))
    } else {
      split_otu4<-rbind(data.frame(tax="___AAA",
                                   value=0,
                                   variable="UP"),
                        split_otu4)
      split_otu4$label<-abs(split_otu4$value)
      color_taxa_use<-c(NA,color_taxa[1:length(unique(microchatLimmaobj$taxon_table[,taxa]))])
      names(color_taxa_use)<-c("___AAA",unique(microchatLimmaobj$taxon_table[,taxa]))
    }

    if ((nrow(split_otu4)-1)>=15) line.x=1.5 else line.x=1

    if (clip_layout)  p3<- ggplot(data = split_otu4,aes(x = substr(tax,start = 4,stop = 5), y = value)) else p3<- ggplot(data = split_otu4,aes(x = tax, y = value))

      p3 + annotate(geom = "rect",inherit.aes=FALSE,check.aes=FALSE,
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                    fill = sigcolor[1],
                    alpha =  sigalpha[2])+
        annotate(geom = "rect",inherit.aes=FALSE,check.aes=FALSE,
                 xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
                 fill = sigcolor[3],
                 alpha =  sigalpha[1])+
      ylim(ifelse(min(split_otu4$value)*1.1>=-3,-3*2,
                  min(split_otu4$value)*1.1*2),max(split_otu4$value)*1.3)+
      ggchicklet::geom_chicklet(data = split_otu4[which(split_otu4$value>0),],
                                aes(fill = tax),
                                width = 0.75,color= adjustcolor(sigcolor[3],0))+
      ggchicklet::geom_chicklet(data = split_otu4[which(split_otu4$value==0),],
                                 fill=NA,
                                  width = 0.75,color= NA)+
      ggchicklet::geom_chicklet(data = split_otu4[which(split_otu4$value<0),],
                                aes(fill = tax),
                                width = 0.75,color= adjustcolor(sigcolor[1],0))+
      scale_fill_manual(values =color_taxa_use,name=taxa)+
      annotate(geom="segment",color="black",lineend="round",size=0.5,
               x = line.x, xend = length(unique(split_otu4$tax))+0.375,
               y = 0, yend = 0)+ylab("OTUs number")+xlab(taxa)+
      theme(panel.grid = element_blank())->p3


    ## add up & down line in the left corner
    p3+annotate(geom="segment",color="black",lineend="round",size=0.5,
                x = line.x, xend = line.x,
                y = 0, yend =  ifelse(max(split_otu4$value)*1.1/2<=3,3*(2-up_coef),
                                      max(split_otu4$value)*1.1/2*(2-up_coef)))+
      annotate(geom="segment",color="black",lineend="round",size=0.5,
               x = line.x, xend = line.x,
               y = ifelse(max(split_otu4$value)*1.1/2<=3,3*up_coef,
                          max(split_otu4$value)*1.1/2*up_coef),
               yend = Inf)+
      geom_point(data=data.frame(x=line.x,y=ifelse(max(split_otu4$value)*1.1/2<=3,3*(2-up_coef),
                                                    max(split_otu4$value)*1.1/2*(2-up_coef))),
                 aes(x,y),fill="black")+
      geom_point(data=data.frame(x=line.x,y=ifelse(max(split_otu4$value)*1.1/2<=3,3*up_coef,
                                                    max(split_otu4$value)*1.1/2*up_coef)),
                 aes(x,y),fill="black")+
      annotate(geom="text",x = line.x, y = ifelse(max(split_otu4$value)*1.1/2<=3,(3*(2-up_coef)+3*up_coef)/2,
                                                   (max(split_otu4$value)*1.1/2*(2-up_coef)+max(split_otu4$value)*1.1/2*up_coef)/2),
               size=2,
               label = "Up",family="serif",fontface="bold")+

      annotate(geom="segment",color="black",lineend="round",size=0.5,
               x = line.x, xend = line.x,
               y = 0, yend = ifelse(min(split_otu4$value)*1.1/2>=-3,-3*(2-dwn_coef),
                                    min(split_otu4$value)*1.1/2*(2-dwn_coef)))+
      annotate(geom="segment",color="black",lineend="round",size=0.5,
               x = line.x, xend = line.x,
               y = ifelse(min(split_otu4$value)*1.1/2>=-3,-3*dwn_coef,
                          min(split_otu4$value)*1.1/2*dwn_coef),
               yend = -Inf)+
      geom_point(data=data.frame(x=line.x,y=ifelse(min(split_otu4$value)*1.1/2>=-3,-3*dwn_coef,
                                                    min(split_otu4$value)*1.1/2*dwn_coef)),
                 aes(x,y),fill="black")+
      geom_point(data=data.frame(x=line.x,y=ifelse(min(split_otu4$value)*1.1/2>=-3,-3*(2-dwn_coef),
                                                    min(split_otu4$value)*1.1/2*(2-dwn_coef))),
                 aes(x,y),fill="black")+
      annotate(geom="text",x = line.x, y = ifelse(min(split_otu4$value)*1.1/2>=-3,((-3)*dwn_coef+(-3)*(2-dwn_coef))/2,
                                                   (min(split_otu4$value)*1.1/2*dwn_coef+min(split_otu4$value)*1.1/2*(2-dwn_coef))/2),
               size=2,
               label = "Down",family="serif",fontface="bold")->p3

    if (clip_layout) {
      p3<-p3+coord_flip()+
        geom_text(aes(label=ifelse(split_otu4$value==0,"",label)),
                  nudge_y=ifelse(split_otu4$value<=0,-1,1),size = 3,color = "black",family = "serif", fontface="bold")+
        theme(aspect.ratio = 1,
              panel.background = element_rect(fill = NA),
              panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5),
              text = element_text(family = 'serif'),
              axis.text.x = element_text(colour='black',size=10,family = "serif"),
              #axis.ticks.y = element_blank(),
              axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif"),
              title = element_text(size = 12),
              axis.text = element_text(face="bold"),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.title.x=element_text(colour='black', size=12,face = "bold",family = "serif"),
              legend.title = element_blank(),
              legend.text=element_text(size=6,face = "bold",colour = "black",
                                       family = "serif"),
              legend.position = legend.position,
              legend.key = element_rect(colour = NA,fill = NA),
              legend.key.width = unit(0.4,"cm"),
              legend.key.height = unit(0.4,"cm"),

              legend.background = element_blank())
    } else {
      p3<-p3+geom_text(aes(label=ifelse(split_otu4$value==0,"",label)),
                       nudge_y=ifelse(split_otu4$value<=0,-1,1),size = 3,color = "black",family = "serif", fontface="bold")+
        theme(aspect.ratio = 1,
              panel.background = element_rect(fill = NA),
              panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5),
              text = element_text(family = 'serif'),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_text(colour='black', size=12,face = "bold",family = "serif"),
              title = element_text(size = 12),
              axis.text = element_text(face="bold"),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.title.y=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
              legend.title = element_blank(),
              legend.text=element_text(size=6,face = "bold",colour = "black",
                                       family = "serif"),
              legend.position = legend.position,
              legend.key = element_rect(colour = NA,fill = NA),
              legend.key.width = unit(0.4,"cm"),
              legend.key.height = unit(0.4,"cm"),
              legend.background = element_blank())
    }
  }
  p3<-p3+labs(title = paste(treat," regulation patterns compared to ",control,sep = ""))+
        theme(title =  element_text(size = 7.5,face="bold"),
              plot.background = element_rect(fill = NA,color = NA))+scale_y_continuous(labels = function(x) abs(x))
  ggsave(paste(export_path,"/DiffTaxaDistributionBarplot_",taxa,"_",paste(control,treat,sep = "-"),"_pvalue(",ptreat$pvalue_thres,")_fc(",sprintf("%0.2f",ptreat$foldchange_thres),").pdf",
               sep = ""),
         units = "cm",
         width = 21/3,
         height = 21*p3$theme$aspect.ratio/3,
         p3)
  return(list(plot=p3,split_otu4=split_otu4,color_taxa_use=color_taxa_use))
}

'plotCombinedDiffTaxaDistri' <- function(mxtreat,
                                         ncol=4,
                                         taxa='Phylum',
                                         sigalpha=cc(0.2,0.2),
                                         color_taxa=colorCustom(50,pal = "gygn"),
                                         dwn_coef=0.4,
                                         up_coef=1.35,
                                         legend.position = c(0.35,0.7),
                                         clip_layout=FALSE,  ##是否旋转下层节点分布图
                                         export_path) {
  mtreat<-mxtreat$data
  xp<-list()
  split_otu_data<-list()
  color_taxa_use_all<-list()
  for (t in 1:length(mtreat)) {
    plotMicrochatDiffTaxaDistri(ptreat=mtreat[[t]],
                                taxa=taxa,
                                sigalpha=sigalpha,
                                color_taxa=color_taxa,
                                dwn_coef=dwn_coef,
                                up_coef=up_coef,
                                legend.position = legend.position,
                                clip_layout=clip_layout,  ##是否旋转下层节点分布图
                                export_path)->ptreat3
    xp[[t]]<-ptreat3$plot
    split_otu_data[[t]]<-ptreat3$split_otu4
    color_taxa_use_all[[t]]<-ptreat3$color_taxa_use
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

  rownum<-length(xp)/ncol
  if (rownum>as.integer(rownum)) nrows<-as.integer(rownum)+1 else nrows<-as.integer(rownum)

  ggsave(paste(export_path,"/data_microbiome/microbial differential analysis/Limma/CombinedDiffTaxaDistriBarplot_",taxa,"_pvalue(",mtreat[[1]]$pvalue_thres,")_fc(",sprintf("%0.2f",mtreat[[1]]$foldchange_thres),").pdf",
               sep = ""),
         units = "cm",
         width = 7*ncol,height = 7*nrows,
         mp)
  return(list(plot=mp,split_otu_data=split_otu_data,color_taxa_use_all=color_taxa_use_all))
}


'plotDiffTaxaDistriLgd' <- function(dtreat,ncol=6,export_path) {

  split_otu_data<-dtreat$split_otu_data
  color_taxa_use_all<-dtreat$color_taxa_use_all

  ### add extra legend
  xxa<-split_otu_data[[1]]$tax
  for (k in 2:length(split_otu_data)) {
    xxa<-union(xxa,split_otu_data[[k]]$tax)
  }
  color_taxa_use<-color_taxa_use_all[[1]][match(xxa,names(color_taxa_use_all[[1]]))]
  color_taxa_use<-color_taxa_use[!is.na(color_taxa_use)]

  color_taxa_use%>%as.data.frame()->extra_p
  extra_p$tax<-rownames(extra_p)
  colnames(extra_p)<-c("color","tax")
  extra_p$tax<-factor(extra_p$tax,levels = unique(extra_p$tax))
  extra_p$value<-0
  ggplot(extra_p, aes(x = tax, y = value)) +
    geom_point(aes(color = tax), size = 0, stroke = 0)+
    scale_color_manual(values = color_taxa_use) +
    guides(color = guide_legend(override.aes = list(
      size = 4, shape = 21,
      fill = extra_p$color),
      direction="horizontal",ncol = ncol,
      label.theme=element_text(size = 6,face="bold",family = "serif"),
      keywidth = unit(0.5,'cm'),
      keyheight = unit(0.5,'cm')))+theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_blank(),
        line= element_blank(),
        rect= element_blank(),
        title= element_blank(),
        legend.position = "top",
        aspect.ratio = 0.001)->lgd_p

  ggsave(paste(export_path,"/data_microbiome/microbial differential analysis/Limma/Legend at ",substr(split_otu_data[[1]]$tax[[3]],start = 1,stop = 3)," levels_for combinedDiffTaxaBarplot.pdf",
               sep = ""),
         units = "cm",
         width = 21,height = 7,
         lgd_p)
  return(lgd_p)
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
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/Limma",sep = "")
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
  #cat("\n","Muti-differential analysis results based on edgeR and limma (",select_comparison,") has been exported to ","/",export_path,"",sep = "","\n")

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
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/Limma",sep = "")
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
  #cat("Muti-volcano plot has been exported to","/",export_path,"",sep = "")
  return(p2)

}

"plotMicrochatDiffCombinedVolcano" <- function(submchat,
                                             ncol=3,
                                             my_comparisons,
                                             xlabname=NULL,
                                             remove.grid=TRUE,
                                             add.barplot=TRUE,
                                             sigcolor=colorCustom(3,pal = "set1"),
                                             color_group=NULL,
                                             sigalpha=c(1,0.3,1),
                                             layout="choice2",
                                             pvalue_thres=0.05,
                                             foldchange_thres=1,
                                             point.size=1,
                                             linetype="dashed",
                                             color_pie.text="black",
                                             color_background=NA,
                                             mytheme=NULL,
                                             export_path) {
if (!is.null(color_group)) color_group<-color_group else sigcolor=sigcolor
  sample.dat<-group_generate(submchat$otu_table)
  g.name<-sample.dat$group%>%unique()
  names(g.name)<-xlabname

  xp<-list()
  mtreat<-list()
  for (i in 1:length(my_comparisons)) {

    #g.name.sel<-my_comparisons[[i]]
    #xlabname.sel<-g.name[match(g.name.sel,g.name)]

    microchatLimmaobj<-calcMicrochatLimma(submchat,
                                          comparison = paste(my_comparisons[[i]][1],
                                                             "-",
                                                             my_comparisons[[i]][2],sep = ""),
                                          file.save=TRUE,
                                          export_path=export_path)

    ### single volcano
    ptreat2<-plotMicrochatDiffVolcano(microchatLimmaobj,
                                      xlabname=xlabname,
                                      remove.grid=remove.grid,
                                      add.barplot=add.barplot,
                                      sigcolor=sigcolor,
                                      color_group=color_group,
                                      sigalpha=sigalpha,
                                      layout=layout,
                                      pvalue_thres=pvalue_thres,
                                      foldchange_thres=foldchange_thres,
                                      point.size=point.size,
                                      linetype=linetype,
                                      color_pie.text=color_pie.text,
                                      color_background=color_background,
                                      mytheme=mytheme,
                                      export_path=export_path)

    xp[[i]]<-ptreat2$plot
    mtreat[[i]]<-ptreat2
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

  rownum<-length(xp)/ncol
  if (rownum>as.integer(rownum)) nrows<-as.integer(rownum)+1 else nrows<-as.integer(rownum)

  ggsave(paste(export_path,"/data_microbiome/microbial differential analysis/Limma/Combined_volcano","_pvalue(",pvalue_thres,")_fc(",sprintf("%0.2f",foldchange_thres),").pdf",
               sep = ""),
         units = "cm",
         width = 7*ncol,height = 7*nrows,
         mp)
return(list(data=mtreat,plot=mp))
}

## lefse----
### lefse====
"transCladeToGeomText" <- function(cladedat,clade_label_level) {
  #4 the 4th try
  tax.ord<-c("p","c","o","f","g","s")
  tax.all.c<-unique(cladedat$node_class)
  tax.all.c<-tax.all.c[-which(tax.all.c=="r" | tax.all.c=="k")]
  tax.ord.re<-tax.ord[ which(tax.all.c %in% tax.ord)]
  cladedat$all.taxa<-NA

  c.check.node<-c()
  for (t in tax.ord.re) {
    pdat<-cladedat[which(cladedat$node_class==t),]
    xdata<-cladedat[which(cladedat$node_class=="p"),]
    if (t=="p") {
      for (sp in 1:nrow(pdat)) {
        cladedat$all.taxa[pdat$node[sp]]<-cladedat$label[pdat$node[sp]]
        c.check.sel<-which(pdat$node[sp] ==cladedat$parent)
        if (length(c.check.sel)!=0) {
          cladedat$all.taxa[c.check.sel]<-paste(cladedat$label[pdat$node[sp]],
                                                cladedat$label[c.check.sel],sep = ";")
          c.check.node<-c(c.check.node,c.check.sel)
        }
      }
    } else {
      for (sp in 1:nrow(pdat)) {
        cladedat$all.taxa[pdat$node[sp]]<-cladedat$all.taxa[pdat$node[sp]]
        c.check.sel<-which(pdat$node[sp] ==cladedat$parent)
        if (length(c.check.sel)!=0) {
          cladedat$all.taxa[c.check.sel]<-paste(cladedat$all.taxa[pdat$node[sp]],
                                                cladedat$label[c.check.sel],sep = ";")
          c.check.node<-c(c.check.node,c.check.sel)
        }
      }
    }
  }
  c.check.node<-c(xdata$node,c.check.node)
  cladedat<-cladedat[-setdiff(cladedat$node,c.check.node),]
  cladedat[,tax.ord]<-str_split_fixed(cladedat$all.taxa,";",length(tax.ord))

  lapply(
    sapply(tax.ord, function(x) {
      cladedat[,x]%>%table()
    }), function(y) {
      tempdat<-as.data.frame(y)
    })->tdata

  ttdat<-tdata[[1]]
  ttdat$Freq<-ttdat$Freq-1
  for (k in 2:length(tdata)) {
    colnames(tdata[[k]])<-colnames(ttdat)
    ttdat<-rbind(ttdat,tdata[[k]])
  }
  colnames(ttdat)[1]<-"label"
  ttdat<-ttdat[-which(ttdat$label==''),]
  rownames(ttdat)<-1:nrow(ttdat)

  cladedat<-merge(cladedat,ttdat,by="label")

  cladedat<-subset(cladedat, select=c(label,x,y,p,c,o,f,g,s,Freq))

  #1
  cladedat$min.y<-0
  cladedat$max.y<-0
  for (t in unique(cladedat$p)) {
    cladedat[which(cladedat$label==t),]$min.y<-min(cladedat[which(cladedat$p==t),]$y)
    cladedat[which(cladedat$label==t),]$max.y<-max(cladedat[which(cladedat$p==t),]$y)
  }

  uq.cdat<-unique(cladedat$c)[ unique(cladedat$c)!='']
  for (t in uq.cdat) {
    cladedat[which(cladedat$label==t),]$min.y<-min(cladedat[which(cladedat$c==t),]$y)
    cladedat[which(cladedat$label==t),]$max.y<-max(cladedat[which(cladedat$c==t),]$y)
  }

  uq.odat<-unique(cladedat$o)[ unique(cladedat$o)!='']
  for (t in uq.odat) {
    cladedat[which(cladedat$label==t),]$min.y<-min(cladedat[which(cladedat$o==t),]$y)
    cladedat[which(cladedat$label==t),]$max.y<-max(cladedat[which(cladedat$o==t),]$y)
  }

  uq.fdat<-unique(cladedat$f)[ unique(cladedat$f)!='']
  for (t in uq.fdat) {
    cladedat[which(cladedat$label==t),]$min.y<-min(cladedat[which(cladedat$f==t),]$y)
    cladedat[which(cladedat$label==t),]$max.y<-max(cladedat[which(cladedat$f==t),]$y)
  }

  uq.gdat<-unique(cladedat$g)[ unique(cladedat$g)!='']
  for (t in uq.gdat) {
    cladedat[which(cladedat$label==t),]$min.y<-min(cladedat[which(cladedat$g==t),]$y)
    cladedat[which(cladedat$label==t),]$max.y<-max(cladedat[which(cladedat$g==t),]$y)
  }
  uq.sdat<-unique(cladedat$s)[ unique(cladedat$s)!='']
  for (t in uq.sdat) {
    cladedat[which(cladedat$label==t),]$min.y<-min(cladedat[which(cladedat$s==t),]$y)
    cladedat[which(cladedat$label==t),]$max.y<-max(cladedat[which(cladedat$s==t),]$y)
  }

  return(cladedat)
}
generate_cladogram_annotation = function(marker_table, tree, color, color_groups, sep = "|") {
  use_marker_table <- marker_table
  feature <- use_marker_table$Taxa
  label <- strsplit(feature, split = sep, fixed = TRUE) %>%
    purrr::map_chr(utils::tail, n =1)
  plot_color <- use_marker_table$group %>%
    as.character
  for(i in seq_along(color_groups)){
    plot_color[plot_color == color_groups[i]] <- color[i]
  }
  annotation <- data.frame(
    node = label,
    color = plot_color,
    enrich_group = use_marker_table$group,
    stringsAsFactors = FALSE
  )
  # filter the feature with bad classification
  annotation %<>% .[label %in% tree$data$label, ]
  annotation
}
get_offset = function(x){(x*0.2+0.2)^2}
get_angle = function(tree, node){
  if (length(node) != 1) {
    stop("The length of `node` must be 1")
  }
  tree_data <- tree$data
  sp <- tidytree::offspring(tree_data, node)$node
  sp2 <- c(sp, node)
  sp.df <- tree_data[match(sp2, tree_data$node),]
  mean(range(sp.df$angle))
}
check_taxa_abund <- function(obj, ...){
  if(is.null(obj$taxa_abund)){
    message("No taxa_abund list found. Calculate it with cal_abund function ...")
    obj$cal_abund(...)
  }
}
check_taxa_level_all = function(taxa_level){
  if(taxa_level == "all"){
    stop("The taxa_level parameter cannot be 'all' for this method! Please provide a taxonomic level, such as 'Genus' !")
  }
}
test_mark = function(dataframe, group, min_num_nonpara = 1, method = NULL){
  d1 <- as.data.frame(t(dataframe))
  taxaname <- colnames(d1)[1]
  d1$group <- group
  colnames(d1)[1] <- "Value"
  formu <- reformulate("group", "Value")
  if(any(table(as.character(group))) < min_num_nonpara){
    list(p_value = NA, med = NA)
  }else{
    if(! is.null(method)){
      method <- match.arg(method, c("wilcox.test", "kruskal.test"))
      if(method == "wilcox.test"){
        res1 <- wilcox.test(formula = formu, data = d1)
      }else{
        res1 <- kruskal.test(formula = formu, data = d1)
      }
    }else{
      if(length(unique(as.character(group))) == 2){
        res1 <- wilcox.test(formula = formu, data = d1)
      }else{
        res1 <- kruskal.test(formula = formu, data = d1)
      }
    }
    if(is.nan(res1$p.value)){
      res1$p.value <- 1
    }
    med <- tapply(d1[,1], group, median) %>% as.data.frame
    colnames(med) <- taxaname
    list(p_value = res1$p.value, med = med)
  }
}
check_taxa_number = function(sel_taxa, p_adjust_method){
  if(sum(sel_taxa) == 0){
    if(p_adjust_method == "none"){
      stop('No significant feature found!')
    }else{
      stop('No significant feature found! To disable p value adjustment, please use p_adjust_method = "none"!')
    }
  }
  if(sum(sel_taxa) == 1){
    if(p_adjust_method == "none"){
      stop('Only one significant feature found! Stop running subsequent process!')
    }else{
      stop('Only one significant feature found! Stop running subsequent process! To disable p value adjustment, please use p_adjust_method = "none"!')
    }
  }
}
generate_p_siglabel <- function(x, nonsig = ""){
  cut(x, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", nonsig))
}
filter_lowabund_feature <- function(abund_table, filter_thres){
  output <- list()
  if(filter_thres > 0){
    mean_abund <- apply(abund_table, 1, mean)
    if(filter_thres > max(mean_abund)){
      stop("Parameter filter_thres is larger than the maximum of mean abundances of features!")
    }
    abund_table <- abund_table[mean_abund >= filter_thres, ]
    filter_features <- mean_abund[mean_abund < filter_thres]
    message("Filter out ", length(filter_features), " features with low abundance ...")
    output[["filter_features"]] <- filter_features
  }
  output[["abund_table"]] <- abund_table
  output
}
dropallfactors <- function(x, unfac2num = FALSE, char2num = FALSE){
  # check x class
  if(!is.data.frame(x)){
    stop("input data must be data.frame class")
  }
  if(unfac2num == T){
    x[] <- lapply(x, function(x) trycharnum(x))
  }else{
    x[] <- lapply(x, function(x) if(is.factor(x)) trycharnum(x) else x)
  }
  if(char2num == T){
    x[] <- lapply(x, function(x) if(is.character(x)) as.factor(x) else x)
    x[] <- lapply(x, function(x) as.numeric(x))
  }
  x
}
trycharnum <- function(x){
  if(suppressWarnings(sum(is.na(as.numeric(as.character(x)))) != sum(is.na(x)))) {
    x <- as.character(x)
  } else {
    x <- as.numeric(as.character(x))
  }
  x
}

calc_diff_lf <- function(datasetx,
                         group = "group",
                         taxa_level = "all",
                         filter_thres = 0,
                         alpha = 0.05,
                         p_adjust_method = "none",
                         lefse_subgroup = NULL,
                         lefse_min_subsam = 10,
                         lefse_norm = 1000000,
                         nresam = 0.6667,
                         boots = 30,
                         remove_unknown = TRUE) {
  datat=list()
  sampleinfo <- datasetx$sample_table
  check_taxa_abund(datasetx)
  if(grepl("all", taxa_level, ignore.case = TRUE)){
    abund_table <- do.call(rbind, unname(datasetx$taxa_abund))
  }else{
    if(! taxa_level %in% names(datasetx$taxa_abund)){
      message("Provided taxa_level: ", taxa_level, " not in tax_table of datasetx; use features in otu_table ...")
      datasetx$add_rownames2taxonomy(use_name = taxa_level)
      suppressMessages(datasetx$cal_abund(rel = TRUE))
    }
    abund_table <- datasetx$taxa_abund[[taxa_level]]
  }

  message(nrow(abund_table), " input features ...")

  all_taxa_input<-nrow(abund_table)
  filter_output <- filter_lowabund_feature(abund_table = abund_table, filter_thres = filter_thres)
  abund_table <- filter_output$abund_table
  filter_features <- filter_output$filter_features

  if(remove_unknown){
    abund_table %<>% {.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.), ignore.case = TRUE), ]}
    message(nrow(abund_table), " features are remained after removing unknown features ...")
    if(nrow(abund_table) == 0){
      stop("No available feature! Please set the parameter: remove_unknown = FALSE")
    }
  }

  abund_table %<>% {. * lefse_norm}
  datat$lefse_norm <- lefse_norm

  group_vec <- sampleinfo[, group] %>% as.factor
  comparisions <- paste0(levels(group_vec), collapse = " - ")
  message("Start Kruskal-Wallis rank sum test for ", group, " ...")
  res_class <- suppressWarnings(lapply(seq_len(nrow(abund_table)), function(x) test_mark(abund_table[x, ], group_vec, method = "kruskal.test")))

  pvalue_raw <- unlist(lapply(res_class, function(x) x$p_value))
  names(pvalue_raw) <- rownames(abund_table)
  pvalue_raw[is.nan(pvalue_raw)] <- 1
  message(sum(pvalue_raw < alpha), " taxa found significant ...")

  pvalue <- p.adjust(pvalue_raw, method = p_adjust_method)
  sel_taxa <- pvalue < alpha
  message("After P value adjustment, ", sum(sel_taxa), " taxa found significant ...")
  check_taxa_number(sel_taxa, p_adjust_method)
  abund_table_sub <- abund_table[sel_taxa, ]
  pvalue_sub <- pvalue[sel_taxa]
  class_taxa_median_sub <- lapply(res_class, function(x) x$med) %>% do.call(cbind, .) %>% .[, sel_taxa]

  all_class_pairs <- combn(unique(as.character(group_vec)), 2)
  # check the difference among subgroups
  if(!is.null(lefse_subgroup)){
    message("Start lefse subgroup biomarkers check for ", lefse_subgroup, " ...")
    all_sub_number <- as.data.table(sampleinfo)[, .N, by = c(group, lefse_subgroup)] %>% as.data.frame %>% .$N
    if(all(all_sub_number < lefse_min_subsam)){
      stop("All sample numbers for subgroups < ", lefse_min_subsam, "! Please consider using small lefse_min_subsam parameter!")
    }
    remove_list_total <- list()
    # for each group paires
    for(i in seq_len(ncol(all_class_pairs))){
      y1 <- all_class_pairs[, i]
      y1_sub_pairs <- expand.grid(
        unique(sampleinfo[sampleinfo[, group] == y1[1], lefse_subgroup]),
        unique(sampleinfo[sampleinfo[, group] == y1[2], lefse_subgroup]),
        stringsAsFactors = FALSE
      ) %>% t
      y1_sub_pairs <- y1_sub_pairs[, unlist(lapply(1:ncol(y1_sub_pairs), function(x){
        ifelse(any(c(sum(sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y1_sub_pairs[1, x]) < lefse_min_subsam,
                     sum(sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y1_sub_pairs[2, x]) < lefse_min_subsam)), FALSE, TRUE)
      })), drop = FALSE]
      if(ncol(y1_sub_pairs) == 0) next
      res_sub_total <- list()
      # check each subgroup pairs under fixed group pair condition
      for(j in 1:ncol(y1_sub_pairs)){
        y2 <- y1_sub_pairs[, j]
        abund_table_sub_y2 <- abund_table_sub[, c(
          rownames(sampleinfo[sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y2[1], ]),
          rownames(sampleinfo[sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y2[2], ])
        )]
        group_vec_sub2 <- c(sampleinfo[sampleinfo[, group] == y1[1] & sampleinfo[, lefse_subgroup] == y2[1], group],
                            sampleinfo[sampleinfo[, group] == y1[2] & sampleinfo[, lefse_subgroup] == y2[2], group])
        res_sub <- lapply(seq_len(nrow(abund_table_sub_y2)), function(x) test_mark(abund_table_sub_y2[x,], group_vec_sub2))
        res_sub_total[[j]] <- res_sub
      }
      raw_median <- class_taxa_median_sub[y1, ] %>%
        {.[1, ] > .[2, ]} %>%
        as.vector
      check_median_sub <- sapply(res_sub_total, function(x) unlist(lapply(x, function(y) {y$med[y1, 1] %>% {.[1] > .[2]}}))) %>% as.data.frame
      check_median_sub[] <- lapply(check_median_sub, function(x) x == raw_median)
      check_p_sub <- sapply(res_sub_total, function(x) unlist(lapply(x, function(y) y$p_value))) %>% as.data.frame
      remove_list <- unlist(lapply(seq_len(nrow(check_median_sub)), function(x){
        if(all(unlist(check_median_sub[x, ]))){
          FALSE
        }else{
          if(any(check_p_sub[x, !unlist(check_median_sub[x, ])] < alpha)){
            TRUE
          }else{
            FALSE
          }
        }
      }))
      remove_list_total[[i]] <- remove_list
    }
    if(!identical(remove_list_total, list())){
      remove_list_total %<>% do.call(cbind, .) %>% apply(., 1, any)
      message("Remove ", sum(remove_list_total), " biomarkers after subgroup check ...")
      abund_table_sub %<>% .[!remove_list_total, ]
      if(nrow(abund_table_sub) == 0){
        stop("No biomarkers remained after subgroup check! stop running!")
      }
      pvalue_sub %<>% .[!remove_list_total]
      class_taxa_median_sub %<>% .[, !remove_list_total]
    }
  }
  res_lda <- list()
  for(num in seq_len(boots)){
    res_lda_pair <- list()
    sample_names_resample <- colnames(abund_table_sub)[base::sample(1:ncol(abund_table_sub), size = ceiling(ncol(abund_table_sub) * nresam))]
    abund_table_sub_resample <- abund_table_sub[, sample_names_resample]
    sampleinfo_resample <- sampleinfo[sample_names_resample, , drop = FALSE]
    # make sure the groups and samples number available
    if(sum(table(as.character(sampleinfo_resample[, group])) > 1) < 2){
      res_lda[[num]] <- NA
      next
    }
    for(i in seq_len(ncol(all_class_pairs))){
      sel_samples <- sampleinfo_resample[, group] %in% all_class_pairs[, i]
      if(length(table(as.character(sampleinfo_resample[sel_samples, group]))) < 2){
        res_lda_pair[[i]] <- NA
        next
      }
      if(min(table(as.character(sampleinfo_resample[sel_samples, group]))) < 2){
        res_lda_pair[[i]] <- NA
        next
      }
      group_vec_lda <- sampleinfo_resample[sel_samples, group] %>% as.character %>% as.factor
      abund_table_sub_lda <- abund_table_sub_resample[, sel_samples]
      abund_table_sub_lda %<>% .[apply(., 1, sd) > 1.0e-10, ]
      if(is.null(lefse_subgroup)){
        abund1 <- cbind.data.frame(t(abund_table_sub_lda), group = group_vec_lda)
      }else{
        subgroup_vec <- sampleinfo_resample[sel_samples, lefse_subgroup] %>% as.character %>% as.factor
        # consider subgroup as an independent variable
        abund1 <- cbind.data.frame(t(abund_table_sub_lda), group = group_vec_lda, lefse_subgroup = subgroup_vec)
      }
      check_res <- tryCatch(mod1 <- MASS::lda(group ~ ., abund1, tol = 1.0e-10), error = function(e) { skip_to_next <- TRUE})
      if(rlang::is_true(check_res)) {
        res_lda_pair[[i]] <- NA
        next
      }else{
        w <- mod1$scaling[,1]
        if(is.null(names(w))){
          names(w) <- rownames(mod1$scaling)
        }
        w_unit <- w/sqrt(sum(w^2))
        w_unit %<>% {.[!grepl("lefse_subgroup", names(.))]}
        ss <- abund1[, !colnames(abund1) %in% c("group", "lefse_subgroup")]
        xy_matrix <- as.matrix(ss)
        LD <- xy_matrix %*% w_unit
        effect_size <- tapply(LD, group_vec_lda, mean) %>% as.vector %>% {.[1] - .[2]} %>% abs
        coeff <- abs(w_unit * effect_size)
        coeff[is.nan(coeff)] <- 0
        names(coeff) %<>% gsub("^`|`$", "", .)
        rres <- mod1$means %>% as.data.frame
        colnames(rres) %<>% gsub("^`|`$", "", .)
        rres <- rres[, rownames(abund_table_sub_lda), drop = FALSE]
        rres1 <- apply(rres, 2, function(x) abs(x[1] - x[2]))
        res_lda_pair[[i]] <- (rres1 + coeff[names(rres1)]) *0.5
      }
    }
    res_lda[[num]] <- res_lda_pair
  }
  res <- sapply(rownames(abund_table_sub), function(k){
    unlist(lapply(seq_len(ncol(all_class_pairs)), function(p){
      unlist(lapply(res_lda, function(x){ x[[p]][k]})) %>% .[!is.na(.)] %>% mean
    })) %>%
      .[!is.na(.)] %>%
      .[!is.nan(.)] %>%
      max
  })
  res <- sapply(res, function(x) {log10(1 + abs(x)) * ifelse(x > 0, 1, -1)})
  output <- cbind.data.frame(group = apply(class_taxa_median_sub, 2, function(x) rownames(class_taxa_median_sub)[which.max(x)]),
                             LDA = res,
                             P.unadj = pvalue_raw[names(pvalue_sub)],
                             P.adj = pvalue_sub)
  output %<>% .[order(.$LDA, decreasing = TRUE), ]
  output <- cbind.data.frame(Comparison = comparisions, Taxa = rownames(output), Method = "LEfSe", output)
  message("Minimum LDA score: ", range(output$LDA)[1], " maximum LDA score: ", range(output$LDA)[2])
  output$Significance <- generate_p_siglabel(output$P.adj, nonsig = "ns")

  # output taxonomic abundance mean and sd for the final res_abund and enrich group finding in metagenomeSeq or ANCOMBC
  res_abund <- reshape2::melt(rownames_to_column(abund_table_sub/lefse_norm, "Taxa"), id.vars = "Taxa")

  # further calculate mean and sd of res_abund with group parameter
  if(!is.null(group)){
    if(group %in% colnames(sampleinfo)){
      colnames(res_abund) <- c("Taxa", "Sample", "Abund")
      res_abund <- suppressWarnings(dplyr::left_join(res_abund, rownames_to_column(sampleinfo), by = c("Sample" = "rowname")))
      res_abund <- summarySE_inter(res_abund, measurevar = "Abund", groupvars = c("Taxa", group))
      colnames(res_abund)[colnames(res_abund) == group] <- "group"
    }
  }
  datat$res_abund <- res_abund
  message('Taxa abundance table is stored in object$res_abund ...')
  if("Factors" %in% colnames(output)){
    output[, "Factors"] %<>% gsub("\\s+$", "", .)
  }
  if(!is.null(output)){
    if("P.adj" %in% colnames(output)){
      if(!"Significance" %in% colnames(output)){
        output$Significance <- generate_p_siglabel(output$P.adj, nonsig = "")
      }
    }
  }
  datat$res_diff <- output
  datat$method <- "lefse"
  datat$taxa_level <- taxa_level
  # save abund_table for the cladogram
  datat$abund_table <- abund_table
  return(datat)
}

calc_lda_filter <- function(datat,lda_thres=1) {
  output <- datat$res_diff
  output<- output[which(output$LDA>lda_thres),]
  datat$res_diff<-output
  message("Selected threshold is ",lda_thres)
  return(datat)
}

'calcTreatAginstControl' <- function(mxtreat,submchat,xlabname) {
  diff_data<-mxtreat$data[[1]]$diff.data

  diff_data<-subset(diff_data,select=c(1))
  diff_data$group<-1
  colnames(diff_data)[ncol(diff_data)]<-xlabname[2]

  for (t in 2:length(mxtreat$data)) {
    diff_data_x<-mxtreat$data[[t]]$diff.data
    diff_data_x<-subset(diff_data_x,select=c(1))
    diff_data_x$group<-1
    colnames(diff_data_x)[ncol(diff_data_x)]<-xlabname[t+1]

    diff_data<-full_join(diff_data,diff_data_x,by="name")
  }

  diff_data[is.na(diff_data)]<-0

  tatxa<-submchat$taxon_table
  tatxa<-tatxa%>%tibble::rownames_to_column(var = "name")
  diff_data<-merge(diff_data,tatxa,by="name")
  diff_data<-subset(diff_data,select=c(1,(ncol(diff_data)-6):ncol(diff_data),2:(ncol(diff_data)-7)))
  ss<-list()
  ss$data_venn3<-diff_data

  diff_data[,c(1,(ncol(diff_data)-length(xlabname)+2):ncol(diff_data))]->shared.data
  shared.data<-shared.data%>%tibble::column_to_rownames(var = "name")

  return(list(venn_data=ss,shared.data=shared.data))
}

cal_abundl = function(datat,
                      select_cols = NULL,
                      rel = TRUE,
                      merge_by = "|",
                      split_group = FALSE,
                      split_by = "&&",
                      split_column = NULL){
  taxa_abund <- list()
  if(is.null(datat$tax_table)){
    stop("No tax_table found! Please check your data!")
  }
  if(nrow(datat$tax_table) != nrow(datat$otu_table)){
    message("The row number of tax_table is not equal to that of otu_table ...")
    message("Automatically applying tidy_datasetx() function to trim the data ...")
    datat$tidy_dataset()
    print(datat)
  }
  if(nrow(datat$sample_table) != ncol(datat$otu_table)){
    message("The sample numbers of sample_table is not equal to that of otu_table ...")
    message("Automatically applying tidy_datasetx() function to trim the data ...")
    datat$tidy_dataset()
    print(datat)
  }
  if(nrow(datat$tax_table) == 0){
    stop("0 row in tax_table! Please check your data!")
  }
  if(is.null(select_cols)){
    select_cols <- seq_len(ncol(datat$tax_table))
  }else{
    if(!is.numeric(select_cols)){
      if(any(! select_cols %in% colnames(datat$tax_table))){
        stop("Part of input names of select_cols are not in the tax_table!")
      }else{
        select_cols <- match(select_cols, colnames(datat$tax_table))
      }
    }
  }
  if(split_group){
    if(is.null(split_column)){
      stop("Spliting rows by one or more columns require split_column parameter! Please set split_column and try again!")
    }
  }
  i=2
  for(i in seq_along(select_cols)){
    taxrank <- colnames(datat$tax_table)[select_cols[i]]
    # assign the columns used for the splitting
    if(!is.null(split_column)){
      if(is.list(split_column)){
        use_split_column <- split_column[[i]]
      }else{
        use_split_column <- split_column
      }
    } else {
      use_split_column <- split_column
    }
#input<-datat
    taxa_abund[[taxrank]] <- microchat::transform_data_proportion(
      datat,
      columns = select_cols[1:i],
      rel = rel,
      merge_by = merge_by,
      split_group = split_group,
      split_by = split_by,
      split_column = use_split_column)
  }
  datat$taxa_abund <- taxa_abund
  message('The result is stored in object$taxa_abund ...')
  return(datat)
}

transform_data_proportion = function(
    input,
    columns,
    rel,
    merge_by = "|",
    split_group = FALSE,
    split_by = "&&",
    split_column){
  library(data.table)
  sampleinfo <- input$sample_table
  abund <- input$otu_table
  tax <- input$tax_table
  tax <- tax[, columns, drop = FALSE]
  # split rows to multiple rows if multiple correspondence
  if(split_group){
    merge_abund <- cbind.data.frame(tax, abund)
    split_merge_abund <- tidyr::separate_rows(merge_abund, all_of(split_column), sep = split_by)
    new_tax <- split_merge_abund[, columns, drop = FALSE]
    new_abund <- split_merge_abund[, (ncol(tax) + 1):(ncol(merge_abund)), drop = FALSE]
    abund1 <- cbind.data.frame(Display = apply(new_tax, 1, paste, collapse = merge_by), new_abund)
  }else{
    abund1 <- cbind.data.frame(Display = apply(tax, 1, paste, collapse = merge_by), abund)
  }
  # first convert table to long format
  abund1 <-abund1 %>%
    as.data.table() %>%
    data.table::melt(id.vars = "Display", value.name = "Abundance", variable.name = "Sample") %>%
    group_by(Display, Sample) %>%
    summarize(sum_abund = sum(Abundance)) %>% data.table() %>%
    #select(-Abundance) %>%
    data.table::setkey(Display, Sample) %>%
    distinct()

  if(rel == T){
    abund1 <- abund1 %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(`res_abund` = sum_abund / sum(sum_abund)) %>%
      dplyr::select(-sum_abund)
  }else{
    colnames(abund1)[colnames(abund1) == "sum_abund"] <- "res_abund"
  }
  abund2 <- abund1 %>%
    as.data.table() %>%
    dcast(Display ~ Sample, value.var = "res_abund") %>%
    as.data.frame() %>%
    `row.names<-`(.[, 1]) %>%
    dplyr::select(-1)
  abund2 <- abund2[order(apply(abund2, 1, mean), decreasing = TRUE), rownames(sampleinfo), drop = FALSE]
  return(abund2)
}

expand_colors <- function(color_values, output_length){
  if(output_length <= length(color_values)){
    total_colors <- color_values[1:output_length]
  }else{
    message("Input colors are not enough to use. Add more colors automatically via color interpolation ...")
    ceiling_cycle_times <- ceiling(output_length/length(color_values))
    total_cycle_times <- lapply(seq_along(color_values), function(x){
      if((ceiling_cycle_times - 1) * length(color_values) + x <= output_length){
        ceiling_cycle_times
      }else{
        ceiling_cycle_times - 1
      }
    }) %>% unlist
    total_color_list <- lapply(seq_along(color_values), function(x){
      colorRampPalette(c(color_values[x], "white"))(total_cycle_times[x] + 1)
    })
    total_colors <- lapply(seq_len(ceiling_cycle_times), function(x){
      unlist(lapply(total_color_list, function(y){
        if((x + 1) <= length(y)){
          y[x]
        }
      }))
    }) %>% unlist
  }
  total_colors
}

generate_microtable_unrel = function(use_datasetx, taxa_level, filter_thres, filter_features){
  use_datasetx$tidy_datasetx()
  suppressMessages(use_datasetx$cal_abund(rel = FALSE))
  use_feature_table <- use_datasetx$taxa_abund[[taxa_level]]
  if(filter_thres > 0){
    use_feature_table %<>% .[! rownames(.) %in% names(filter_features), ]
  }
  message("Available feature number: ", nrow(use_feature_table))
  newdata <- microtable$new(otu_table = use_feature_table, sample_table = use_datasetx$sample_table)
  newdata$tidy_datasetx()
  newdata
}

##plot
plot_backgroud_tree = function(abund_table, use_taxa_num = NULL, filter_taxa = NULL, sep = "|",branch_size=0.2,layout = 'circular',
                               same.branch.size=TRUE,color_TaxaOrClass=TRUE,
                               branch.size=c(4,3.5,3,2.5,2,1.5,1,0.5),
                               color_branch=colorCustom(8,pal = "pkgn")){
  # filter the taxa with unidentified classification or with space, in case of the unexpected error in the following operations
  abund_table %<>% {.[!grepl("\\|.__\\|", rownames(.)), ]} %>%
    {.[!grepl("\\s", rownames(.)), ]} %>%
    # also filter uncleared classification to make it in line with the lefse above
    {.[!grepl("Incertae_sedis|unculture", rownames(.), ignore.case = TRUE), ]}
  if(nrow(abund_table) <= 2){
    stop("After filtering out non-standard taxonomy information, the abundance table only has ", nrow(abund_table), " feature(s)! ",
         "Is there an issue with the taxonomy table? ",
         "Please first use the tidy_taxonomy function to process the taxonomy information table before constructing the microtable object.")
  }
  if(!is.null(use_taxa_num)){
    if(use_taxa_num < nrow(abund_table)){
      message("Select ", use_taxa_num, " most abundant taxa as the background cladogram ...")
      abund_table %<>% .[names(sort(apply(., 1, mean), decreasing = TRUE)[1:use_taxa_num]), ]
    }else{
      message("Provided use_taxa_num: ", use_taxa_num, " >= ", " total effective taxa number. Skip the selection ...")
    }
  }
  if(!is.null(filter_taxa)){
    abund_table %<>% .[apply(., 1, mean) > (datat$lefse_norm * filter_taxa), ]
  }
  abund_table %<>% .[sort(rownames(.)), ]
  tree_table <- data.frame(taxa = row.names(abund_table), abd = rowMeans(abund_table), stringsAsFactors = FALSE) %>%
    dplyr::mutate(taxa = paste("r__Root", .data$taxa, sep = sep), abd = .data$abd/max(.data$abd)*100)
  taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
  nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
  # check whether some nodes duplicated from bad classification
  if(any(duplicated(nodes))){
    del <- nodes %>% .[duplicated(.)] %>% unique
    for(i in del){
      tree_table %<>% .[!grepl(paste0("\\|", i, "($|\\|)"), .$taxa), ]
    }
    taxa_split <- strsplit(tree_table$taxa, split = sep, fixed = TRUE)
    nodes <- purrr::map_chr(taxa_split, utils::tail, n = 1)
  }
  # add root node
  nodes %<>% c("r__Root", .)
  # levels used for extending clade label
  label_levels <- purrr::map_chr(nodes, ~ gsub("__.*$", "", .x)) %>%
    factor(levels = rev(unlist(lapply(taxa_split, function(x) gsub("(.)__.*", "\\1", x))) %>% .[!duplicated(.)]))

  # root must be a parent node
  nodes_parent <- purrr::map_chr(taxa_split, ~ .x[length(.x) - 1]) %>% c("root", .)

  # add index for nodes
  is_tip <- !nodes %in% nodes_parent
  index <- vector("integer", length(is_tip))
  index[is_tip] <- 1:sum(is_tip)
  index[!is_tip] <- (sum(is_tip) + 1):length(is_tip)

  edges <- cbind(parent = index[match(nodes_parent, nodes)], child = index)
  edges <- edges[!is.na(edges[, 1]), ]
  # not label the tips
  node_label <- nodes[!is_tip]
  phylo <- structure(list(
    edge = edges,
    node.label = node_label,
    tip.label = nodes[is_tip],
    edge.length = rep(1, nrow(edges)),
    Nnode = length(node_label)
  ), class = "phylo")
  mapping <- data.frame(
    node = index,
    abd = c(100, tree_table$abd),
    node_label = nodes,
    fulltaxa = c(paste0("r__Root","|k__Bacteria|p__total"),tree_table$taxa),
    stringsAsFactors = FALSE)
  mapping$fulltaxa[2]<-mapping$fulltaxa[1]
  mapping$node_class <- label_levels
  mapping$taxa<-0
  for (t in 1:nrow(mapping)) {
    mapping$taxa[t]<- strsplit(mapping$fulltaxa,sep,2)[[t]][3]
  }

  if (!color_TaxaOrClass) {
    colors_class<-color_branch
    names(colors_class)<-unique(mapping$node_class)%>%sort()

    tax.class<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    if (same.branch.size) {
      tree <- tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
      tree <- ggtree::ggtree(tree, size = branch_size, layout = layout,aes(color=node_class))+
        scale_color_manual(values = colors_class,guide=FALSE)
    } else {
      b.size<-branch.size
      node_clas<-substr(tolower(tax.class),start = 1,stop = 1)
      node_clas<-c("r",node_clas)
      names(b.size)<-node_clas
      mapping$bsize <- b.size[match(mapping$node_class,names(b.size))]
      tree <- tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
      tree <- ggtree::ggtree(tree,aes(size=I(bsize),color=node_class), layout = layout)+
        scale_color_manual(values = colors_class,guide=FALSE)
    }
  } else {
    colors_class<-color_branch
    if (length(colors_class)<length(unique(mapping$taxa))) {
      message("Totally displayed ",length(unique(mapping$taxa)), " phyla")
      stop("Upon set the `color_TaxaOrClass` argument to `True`, please provide more colors in the `color_branch` argument !!!")
    }
    names(colors_class)<-unique(mapping$taxa)

    tax.class<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    if (same.branch.size) {
      tree <- tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
      tree <- ggtree::ggtree(tree, size = branch_size, layout = layout,
                             aes(color=taxa))+
        scale_color_manual(values = colors_class,guide=FALSE)
    } else {
      b.size<-branch.size
      node_clas<-substr(tolower(tax.class),start = 1,stop = 1)
      node_clas<-c("r",node_clas)
      names(b.size)<-node_clas
      mapping$bsize <- b.size[match(mapping$node_class,names(b.size))]
      tree <- tidytree::treedata(phylo = phylo, data = tibble::as_tibble(mapping))
      tree <- ggtree::ggtree(tree,aes(size=I(bsize),color=taxa), layout = layout)+
        scale_color_manual(values = colors_class,guide=FALSE)
    }
  }
 return(tree)
}

plot_diff_cladogram = function(datat,
                               color=colorCustom(3,pal = "na5"),
                               color_TaxaOrClass=TRUE,
                               group_order = NULL,
                               use_taxa_num = 100,
                               use_feature_num = 30,
                               select_show_labels = NULL,
                               filter_taxa = NULL,
                               clade_label_level = 4,
                               only_select_show = FALSE,
                               sep = "|",
                               branch_size = 0.4,
                               layout = 'circular',
                               alpha = 0.4,
                               pt.stroke=0.4,
                               clade_label_size = 1,
                               clade_label_size_add = 3,
                               clade_label_size_log = exp(2),
                               node_size_scale = 1,
                               node_size_offset = 1,
                               annotation_shape = 22,
                               same.point.shape= TRUE,
                               point.shape = 21:25,
                               annotation_shape_size = 5,
                               same.branch.size=FALSE,
                               branch.size=c(1.4,1.2,1.0,0.8,0.6,0.4,0.2,0.1),
                               color_branch=colorCustom(8,pal = "pkgn"),
                               export_path){
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/LEfSe",sep = "")
  dir.create(export_path, recursive = TRUE)
  tax.class<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  message("Argument 'clade_label_level' you have chosen is ",clade_label_level)
  tax.class.sel<-c()
  for (t in 1:length(tax.class[1:(8-clade_label_level)])) {
    tax.class.sel<-c(tax.class.sel,tax.class[t])
  }
  tax.class.sel<-tax.class.sel[-1]
  tax.class.sel.x<-tax.class.sel[1]
  for (t in 2:length(tax.class.sel)) {
    tax.class.sel.x<-paste(tax.class.sel.x,tax.class.sel[t],sep = ", ")
  }
  message("The full name of the taxon whose classification level is higher than (contains) ", tax.class[8-clade_label_level]," will be displayed, namely ")
  cat(tax.class.sel.x,"\n")
  # developed based on microbiomeMarker
  abund_table <- datat$abund_table
  marker_table <- datat$res_diff %>% dropallfactors
  method <- datat$method

  if(! method %in% c("lefse", "rf")){
    stop("This function currently can only be used for method = 'lefse' or 'rf' !")
  }
  if(datat$taxa_level != "all"){
    stop("This function is available only when taxa_level = 'all' !")
  }
  if(!is.null(use_feature_num)){
    if(use_feature_num > nrow(marker_table)){
      message("Input use_feature_num ", use_feature_num, " larger than available features number ", nrow(marker_table), " ! ",
              "Use ", nrow(marker_table), " instead of it ...")
    }else{
      message("Select ", use_feature_num, " significant features ...")
      marker_table %<>% .[1:use_feature_num, ]
    }
  }
  if(only_select_show == T){
    marker_table %<>% .[.$Taxa %in% select_show_labels, ]
  }
  # color legend order settings
  if(is.null(group_order)){
    if(! is.null(datat$group_order)){
      color_groups <- datat$group_order
    }else{
      color_groups<-unique(substr(colnames(datat$abund_table),start = 1,stop = 2))
      color_groups_ <- marker_table$group %>% as.character %>% as.factor %>% levels
      color_groups<-color_groups[which(color_groups %in% color_groups_)]
    }
  }else{
    color_groups <- group_order
  }

  # filter redundant groups
  if(! all(color_groups %in% unique(marker_table$group))){
    tmp_message <- color_groups[! color_groups %in% unique(marker_table$group)]
    message("Part of groups in group_order, ", paste(tmp_message, collapse = " "), ", not found in the feature table ...")
    color_groups %<>% .[. %in% unique(marker_table$group)]
  }
  if(! all(unique(marker_table$group) %in% color_groups)){
    tmp_message <- unique(marker_table$group)[! unique(marker_table$group) %in% color_groups]
    message("Part of groups in the feature table, ", paste(tmp_message, collapse = " "), ", not found in the group_order ...")
    marker_table %<>% .[.$group %in% color_groups, ]
  }
  # get the color palette
  color_groups_all<-unique(substr(colnames(datat$abund_table),start = 1,stop = 2))
  if(length(color) < length(unique(marker_table$group))){
    stop("Please provide enough color palette! There are ", length(unique(marker_table$group)),
         " groups, but only ", length(color), " colors provideed in color parameter!")
  }else{
    color <- color[match(color_groups,color_groups_all)]
  }

  tree <- plot_backgroud_tree(abund_table = abund_table, use_taxa_num = use_taxa_num,
                              filter_taxa = filter_taxa, sep = sep,color_TaxaOrClass=color_TaxaOrClass,
                              branch_size=branch_size,layout = layout,
                              same.branch.size=same.branch.size,
                              branch.size=branch.size,
                              color_branch=color_branch)
  cladedat<-tree$data%>%as.tibble()
  cladedat<-transCladeToGeomText(cladedat,clade_label_level)
  # generate annotation
  annotation <- generate_cladogram_annotation(marker_table, tree = tree, color = color, color_groups = color_groups, sep = sep)
  # check again filtered groups
  if(!is.null(use_taxa_num)){
    if(! all(color_groups %in% unique(annotation$enrich_group))){
      tmp_message <- color_groups[! color_groups %in% unique(annotation$enrich_group)]
      message("Biomarkers in group(s), ", paste(tmp_message, collapse = " "), ", not found in the background tree!",
              " Please try to enlarge parameter use_taxa_num from ", use_taxa_num, " to a larger one ...")
      color_groups %<>% .[. %in% unique(annotation$enrich_group)]
    }
  }

  annotation_info <- dplyr::left_join(annotation, tree$data, by = c("node" = "label")) %>%
    dplyr::mutate(label = .data$node, id = .data$node.y, level = as.numeric(.data$node_class))
  hilight_para <- dplyr::transmute(
    annotation_info,
    node = .data$id,
    fill = .data$color,
    alpha = alpha,
    extend = get_offset(.data$level)
  )
  hilights_g <- purrr::pmap(hilight_para, ggtree::geom_hilight)
  tree <- purrr::reduce(hilights_g, `+`, .init = tree)

  # hilight legend
  hilights_df <- dplyr::distinct(annotation_info, .data$enrich_group, .data$color)
  hilights_df$x <- 0
  hilights_df$y <- 1
  # resort the table used for the legend color and text
  hilights_df %<>% `row.names<-`(.$enrich_group) %>% .[color_groups, ]
  # make sure the right order in legend
  hilights_df$enrich_group %<>% factor(., levels = color_groups)

  # add legend
  tree <- tree +
    geom_rect(aes(xmin = x, xmax = x, ymax = y, ymin = y, fill = enrich_group), data = hilights_df, inherit.aes = FALSE) +
    guides(fill = guide_legend(title = NULL, order = 1, override.aes = list(fill = hilights_df$color)))

  # set nodes color and size
  nodes_colors <- rep("white", nrow(tree$data))
  nodes_colors[annotation_info$id] <- annotation_info$color
  node_size <- node_size_scale*log(tree$data$abd) + node_size_offset
  tree$data$node_size <- node_size

  if (!same.point.shape) {
    node_clas<-substr(tolower(tax.class),start = 1,stop = 1)
    node_clas<-c("r",node_clas)
    if (length(point.shape)<length(node_clas)) point.shape<-rep(point.shape,7)[1:length(node_clas)]
    names(point.shape)<-node_clas
    tree$data$shape <- point.shape[match(tree$data$node_class,names(point.shape))]
    tree <- tree +
      ggnewscale::new_scale("size")+
      ggtree::geom_point2(
        aes(size = I(node_size)/1.5,shape = I(shape)),
        stroke=pt.stroke,
        fill = nodes_colors)
  } else {
    tree <- tree +
      ggnewscale::new_scale("size")+
      ggtree::geom_point2(
        aes(size = I(node_size)/1.5),
        stroke=pt.stroke,
        fill = nodes_colors, shape = annotation_shape)
  }

  ## add clade labels
  clade_label <- dplyr::transmute(
    annotation_info,
    node = .data$id,
    offset = get_offset(.data$level)-0.4,
    offset.text = 0,
    angle = purrr::map_dbl(.data$id, get_angle, tree = tree),
    label = .data$label,
    fontsize = clade_label_size + log(.data$level + clade_label_size_add, base = clade_label_size_log),
    barsize = 0,
    extend = 0.2,
    hjust = 0.5,
    level = .data$level
  ) %>% dplyr::arrange(desc(.data$level))

  clade_label$offset.text <- unlist(lapply(seq_len(nrow(clade_label)), function(x){
    if(clade_label$angle[x] < 180){ 0.2 }else{ 0 }}))
  clade_label$angle <- unlist(lapply(clade_label$angle, function(x){
    if(x < 180){ x - 90 }else{ x + 90 }}))
  clade_label_new <- clade_label

  # add letters label to replace long taxonomic label
  if(is.null(select_show_labels)){
    # outer circle --> larger level; label smaller levels, i.e. finer taxonomy
    ind <- clade_label$level < clade_label_level
  }else{
    ind <- ! clade_label$label %in% select_show_labels
  }
  ind_num <- sum(ind)

  if(ind_num > 0){
    if(ind_num < 27){
      use_letters <- letters
    }else{
      if(ind_num < 326){
        use_letters <- apply(combn(letters, 2), 2, function(x){paste0(x, collapse = "")})
      }else{
        stop("Too much features to be labelled with letters, consider to use use_feature_num parameter to reduce the number!")
      }
    }
    clade_label_new$label_legend <- clade_label_new$label_show <- clade_label_new$label_raw <- clade_label_new$label
    clade_label_new$label_show[ind] <- use_letters[1:ind_num]
    clade_label_new$label_legend[ind] <- paste0(clade_label_new$label_show[ind], ": ", clade_label_new$label[ind])
    clade_label_new$label <- clade_label_new$label_show
    # delete redundant columns to avoid warnings
    clade_label <- clade_label_new %>% .[, which(! colnames(.) %in% c("label_raw", "label_show", "label_legend", "level"))]
  }

  clade_label_g <- purrr::pmap(clade_label, ggtree::geom_cladelabel)
  #tree <- purrr::reduce(clade_label_g, `+`, .init = tree)

  clade_label_x<-merge(clade_label,tree$data,by="node")
  clade_label_x<-subset(clade_label_x,select=c(2,3,4,20,6,7,8,9,19,
                                               5,12,1,10,15,17,18,21 ))

  cladedat_new<-subset(cladedat,select=c(label,Freq,min.y,max.y))
  colnames(cladedat_new)[1]<-"label.y"
  clade_label_x<-merge(clade_label_x,cladedat_new,by="label.y")
  clade_label_x$ypos<-(clade_label_x$min.y+clade_label_x$max.y)/2

  tree$data->sac
  for (t in 1:nrow(clade_label_x)) {
    clade_label_x$angle[t]<-sac[which(sac$node_label==clade_label_x$label.y[t]),]$angle
  }


  if (layout %in% c("circular","radial","fan",'inward_circular')) {
    if (layout %in% c("radial"))     tree<-tree+
        geom_text(data=clade_label_x,inherit.aes = FALSE,
                         aes(x= ifelse(x<=(8-clade_label_level),
                                       microchat::normalize(9-x)$x+7,
                                       x+0.2),
                             y=ypos,label=label.x,
                             angle=angle.x,
                             size=fontsize),
                         vjust=0,family="serif")

    if (layout %in% c("circular","inward_circular","fan"))   tree<-tree+
        geom_text(data=clade_label_x,inherit.aes = FALSE,
                   aes(x= ifelse(x<=(8-clade_label_level),
                                 microchat::normalize(9-x)$x+7,
                                 x+0.2),
                       y=ypos,label=label.x,
                       angle=angleAdj(angle),
                       size=fontsize),
                   vjust=0,family="serif")
  } else {
    tree<-tree+
      geom_text(data=clade_label_x,inherit.aes = FALSE,
                aes(x= ifelse(x<=(8-clade_label_level),
                              microchat::normalize(9-x)$x+7,
                              x+0.2),
                    y=ypos,label=label.x,
                    #angle=angle.y,
                    size=fontsize),angle=0,
                vjust=0,family="serif")
  }

  # if letters are used, add guide labels
  if(ind_num > 0){
    guide_label <- clade_label_new[ind, ] %>%
      dplyr::mutate(color = annotation_info$color[match(.data$label_raw, annotation_info$label)])

    if (!same.point.shape) {

      guide_label$shape<- point.shape[match(substr(guide_label$label_legend,start = 4,stop = 4),
                                            names(point.shape))]
      tree <- tree
    } else {
      tree <- tree +
        geom_point(data = guide_label, inherit.aes = FALSE, aes_(x = 0, y = 0, shape = ~label_legend), size = 0, stroke = 0) +
        scale_shape_manual(values = rep(annotation_shape, nrow(guide_label))) +
        guides(shape = guide_legend(override.aes = list(
          size = annotation_shape_size, shape = annotation_shape,
          fill = guide_label$color),
          label.theme=element_text(size = 4,face="bold",family = "serif"),
          keywidth = unit(0.1,'cm'),
          keyheight = unit(0.1,'cm')))
    }
  }
  tree <- tree + theme(legend.position = "right", legend.title = element_blank(),
                       text = element_text(family = "serif"),
                       legend.key.size = unit(0.4,'cm'),
                       legend.text =  element_text(size = 4,face="bold"))

  ggsave(paste(export_path,"/microchat_lefse cladogram",".pdf",sep = ""),
         units = "cm",
         width =21/3*2,
         height =  21/3*2,
         tree)

  return(tree)
}

"angleAdj" <- function(angle) {
  ifelse(angle<45,
         angle-105,
         ifelse(angle<90,
                angle-120,
                ifelse(angle<135,
                       angle-135,
                       ifelse(angle<180,
                              angle-150,
                              ifelse(angle<225,
                                     angle-135,
                                     ifelse(angle<270,
                                            angle-120,
                                            ifelse(angle<315,
                                                   angle-105,
                                                   ifelse(angle<360,
                                                          angle-90,
                                                          angle-75))
                                     )
                              )
                       )
                )
         )
  )
}
plot_diff_abund = function(datat,
                           use_number = 1:20,
                           color_values = RColorBrewer::brewer.pal(8, "Dark2"),
                           select_group = NULL,
                           select_taxa = NULL,
                           simplify_names = TRUE,
                           keep_prefix = TRUE,
                           group_order = NULL,
                           barwidth = 0.9,
                           use_se = TRUE,
                           add_sig = FALSE,
                           add_sig_label = "Significance",
                           add_sig_label_color = "black",
                           add_sig_tip_length = 0.01,
                           y_start = 1.01,
                           y_increase = 0.05,
                           text_y_size = 10,
                           coord_flip = TRUE,
                           xtext_angle = 45,
                           ...){
  abund_data <- datat$res_abund
  method <- datat$method
  diff_data <- datat$res_diff
  if(grepl("ancombc2", method)){
    stop("The function can not be applied to ancombc2!")
  }
  if(grepl("formula", method)){
    stop("The function can not be applied to multi-factor analysis!")
  }
  # first determine how to select compared groups
  if(!is.null(select_group)){
    if(length(select_group) > 1){
      stop("The select_group parameter should only have one element! Please check the input!")
    }
    if(! select_group %in% diff_data$Comparison){
      stop("The select_group parameter must be one of elements of object$res_diff$Comparison!")
    }
    diff_data %<>% .[.$Comparison %in% select_group, ]
    select_group_split <- strsplit(select_group, split = " - ") %>% unlist
    abund_data %<>% .[.$group %in% select_group_split, ]
  }
  # sort according to different columns
  if(method == "metastat"){
    message('Reorder taxa according to qvalue in res_diff from low to high ...')
    diff_data %<>% .[order(.$qvalue, decreasing = FALSE), ]
    # diff_data %<>% .[.$qvalue < 0.05, ]
  }else{
    # lefse and rf are ordered
    if(! method %in% c("lefse", "rf", "anova")){
      if("P.adj" %in% colnames(diff_data)){
        message('Reorder taxa according to P.adj in res_diff from low to high ...')
        diff_data %<>% .[order(.$P.adj, decreasing = FALSE), ]
      }
    }
  }
  if(nrow(diff_data) == 0){
    stop("No significant taxa can be used to plot the abundance!")
  }
  if(simplify_names == T){
    diff_data$Taxa %<>% gsub(".*\\|", "", .)
    abund_data$Taxa %<>% gsub(".*\\|", "", .)
  }
  if(keep_prefix == F){
    diff_data$Taxa %<>% gsub(".__", "", .)
    abund_data$Taxa %<>% gsub(".__", "", .)
  }
  if(is.null(select_taxa)){
    if(length(use_number) > length(unique(as.character(diff_data$Taxa)))){
      message("The length of use_number is larger than taxa number in object$diff_data. Use all taxa ...")
      use_number <- 1:length(unique(as.character(diff_data$Taxa)))
    }
    diff_data %<>% .[.$Taxa %in% unique(as.character(diff_data$Taxa))[use_number], ]
    diff_data$Taxa %<>% factor(., levels = rev(unique(as.character(.))))
  }else{
    message('Use provided select_taxa to filter and reorder taxa ...')
    diff_data %<>% .[.$Taxa %in% select_taxa, ]
    if(nrow(diff_data) == 0){
      stop("No significant taxa can be used to plot the abundance!")
    }
    diff_data$Taxa %<>% factor(., levels = rev(select_taxa))
  }
  abund_data %<>% .[.$Taxa %in% levels(diff_data$Taxa), ]
  abund_data$Taxa %<>% factor(., levels = levels(diff_data$Taxa))
  if(is.null(group_order)){
    if((!is.null(datat$group_order)) & (length(unique(abund_data$group)) == length(datat$group_order))){
      abund_data$group %<>% factor(., levels = rev(datat$group_order))
    }else{
      abund_data$group %<>% as.character %>% as.factor
    }
  }else{
    abund_data$group %<>% factor(., levels = rev(group_order))
  }
  if(length(color_values) < length(levels(abund_data$group))){
    stop("Please provide color_values parameter with more colors!")
  }else{
    color_values %<>% .[1:length(levels(abund_data$group))] %>% rev
  }
  if(!coord_flip){
    abund_data$group %<>% factor(., levels = rev(levels(.)))
    abund_data$Taxa %<>% factor(., levels = rev(levels(.)))
    diff_data$Taxa %<>% factor(., levels = rev(levels(.)))
    color_values %<>% rev
  }

  p <- ggplot(abund_data, aes(x = Taxa, y = Mean, color = group, fill = group)) +
    theme_bw() +
    geom_bar(stat="identity", position = position_dodge(), width = barwidth)
  if(use_se == T){
    p <- p + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.45, position=position_dodge(barwidth), color = "black")
  }else{
    p <- p + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.45, position=position_dodge(barwidth), color = "black")
  }

  if(add_sig){
    # assign labels by factor orders
    x_axis_order <- levels(abund_data$group)
    if(! add_sig_label %in% colnames(diff_data)){
      stop("The add_sig_label parameter must be one of colnames of object$res_diff!")
    }
    if(is.factor(diff_data[, add_sig_label])){
      diff_data[, add_sig_label] %<>% as.character
    }else{
      if(is.numeric(diff_data[, add_sig_label])){
        diff_data[, add_sig_label] %<>% round(., 4)
      }
    }
    if(use_se){
      y_start_use <- max((abund_data$Mean + abund_data$SE)) * y_start
    }else{
      y_start_use <- max((abund_data$Mean + abund_data$SD)) * y_start
    }
    all_taxa <- levels(abund_data$Taxa)

    # for groups > 3, show the global comparision
    if((length(levels(abund_data$group)) > 2 & method %in% c("lefse", "rf", "KW", "ALDEx2_kw")) |
       (length(unlist(gregexpr(" - ", diff_data$Comparison[1]))) > 1 & method == "ancombc2")){
      add_letter_text <- diff_data[match(all_taxa, diff_data$Taxa), add_sig_label]
      textdf <- data.frame(
        x = all_taxa,
        y = y_start_use,
        add = add_letter_text,
        stringsAsFactors = FALSE
      )
      p <- p + geom_text(aes(x = x, y = y, label = add), data = textdf, inherit.aes = FALSE)
    }else{
      if(! "Letter" %in% colnames(diff_data)){
        if(any(grepl("\\s-\\s", x_axis_order))){
          stop("The group names have ' - ' characters, which can hinder the group recognition and mapping in the plot! Please rename groups and rerun!")
        }
        annotations <- c()
        x_min <- c()
        x_max <- c()
        y_position <- c()

        start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
        increase_bar_mid <- barwidth/length(x_axis_order)

        for(j in all_taxa){
          select_use_diff_data <- diff_data %>% dropallfactors %>% .[.$Taxa == j, ]
          for(i in seq_len(nrow(select_use_diff_data))){
            # first determine the bar range
            mid_num <- match(j, all_taxa) - 1
            annotations %<>% c(., select_use_diff_data[i, add_sig_label])
            x_min %<>% c(., mid_num +
                           (start_bar_mid + (match(gsub("(.*)\\s-\\s(.*)", "\\1", select_use_diff_data[i, "Comparison"]), x_axis_order) - 1) * increase_bar_mid))
            x_max %<>% c(., mid_num +
                           (start_bar_mid + (match(gsub("(.*)\\s-\\s(.*)", "\\2", select_use_diff_data[i, "Comparison"]), x_axis_order) - 1) * increase_bar_mid))
            y_position %<>% c(., y_start_use * (1 + i * y_increase))
          }
        }
        p <- p + ggsignif::geom_signif(
          annotations = annotations,
          y_position = y_position,
          xmin = x_min,
          xmax = x_max,
          color = add_sig_label_color,
          tip_length = add_sig_tip_length,
          ...
        )
      }else{
        x_mid <- c()
        annotations <- c()
        y_position <- c()

        start_bar_mid <- 1 - (barwidth/2 - barwidth/(length(x_axis_order) * 2))
        increase_bar_mid <- barwidth/length(x_axis_order)

        for(j in all_taxa){
          select_use_diff_data <- diff_data %>% dropallfactors %>% .[.$Taxa == j, ]
          for(i in seq_len(nrow(select_use_diff_data))){
            mid_num <- match(j, all_taxa) - 1
            annotations %<>% c(., select_use_diff_data[i, add_sig_label])
            x_mid %<>% c(., mid_num + (start_bar_mid + (match(select_use_diff_data[i, "group"], x_axis_order) - 1) * increase_bar_mid))
            abund_data_select <- abund_data[abund_data$group == select_use_diff_data[i, "group"] & abund_data$Taxa == j, ]
            if(use_se){
              y_position %<>% c(., y_start * (abund_data_select$Mean + abund_data_select$SE) + y_increase * max(abund_data$Mean))
            }else{
              y_position %<>% c(., y_start * (abund_data_select$Mean + abund_data_select$SD) + y_increase * max(abund_data$Mean))
            }
          }
        }
        textdf <- data.frame(
          x = x_mid,
          y = y_position,
          add = annotations,
          stringsAsFactors = FALSE
        )
        p <- p + geom_text(aes(x = x, y = y, label = add), data = textdf, inherit.aes = FALSE)
      }
    }
  }

  p <- p + scale_color_manual(values = color_values) +
    scale_fill_manual(values = color_values) +
    ylab("Relative abundance") +
    theme(legend.position = "right") +
    theme(panel.border = element_blank(), panel.background=element_rect(fill="white")) +
    theme(axis.title = element_text(size = 17))

  if(coord_flip){
    p <- p + coord_flip() + guides(fill = guide_legend(reverse = TRUE, ncol = 1), color = "none") +
      theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) +
      theme(axis.title.y = element_blank(), axis.text.y = element_text(size = text_y_size, color = "black"))
  }else{
    p <- p + guides(fill = guide_legend(reverse = FALSE, ncol = 1), color = "none")
    if(xtext_angle != 0){
      p <- p + theme(axis.text.x = element_text(angle = xtext_angle, colour = "black", vjust = 1, hjust = 1, size = text_y_size))
    }else{
      p <- p + theme(axis.text.x = element_text(angle = xtext_angle, colour = "black", size = text_y_size))
    }
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = text_y_size, color = "black")) +
      theme(plot.margin = unit(c(.1, .1, .1, 1), "cm"))
  }
  p
}

plot_diff_bar = function(datat,
                         color_values = RColorBrewer::brewer.pal(8, "Dark2"),
                         color_group_map = FALSE,
                         use_number = 1:10,
                         threshold = NULL,
                         select_group = NULL,
                         keep_full_name = FALSE,
                         keep_prefix = TRUE,
                         group_order = NULL,
                         axis_text_y = 12,
                         coord_flip = TRUE,
                         xtext_angle = 45,
                         xtext_size = 10,
                         heatmap_cell = "P.unadj",
                         heatmap_sig = "Significance",
                         heatmap_x = "Factors",
                         heatmap_y = "Taxa",
                         heatmap_lab_fill = "P value",
                         export_path,
                         ...){
  export_path<-paste(export_path,"/data_microbiome/microbial differential analysis/LEfSe",sep = "")
  dir.create(export_path, recursive = TRUE)
  use_data <- datat$res_diff
  method <- datat$method
  if(is.null(method)){
    if(is.null(use_data)){
      stop("The res_diff should be provided when method is NULL!")
    }else{
      if(!"Value" %in% colnames(use_data)){
        stop("No Value column found in the object$res_diff!")
      }
    }
  }
  if(keep_full_name == F){
    if(any(grepl("\\..__", use_data$Taxa))){
      use_data$Taxa %<>% gsub(".*(.__.*?$)", "\\1", .)
    }else{
      use_data$Taxa %<>% gsub(".*\\|", "", .)
    }
  }
  if(keep_prefix == F){
    use_data$Taxa %<>% gsub(".__", "", .)
  }

    if(method == "lefse"){
      colnames(use_data)[colnames(use_data) == "LDA"] <- "Value"
      ylab_title <- "LDA score"
    }else{
      if(method == "rf"){
        colnames(use_data)[colnames(use_data) == "MeanDecreaseGini"] <- "Value"
        ylab_title <- "MeanDecreaseGini"
      }else{
        if(method == "metastat"){
          use_data %<>% .[.$qvalue < 0.05, ]
          use_data$Value <- 1 - use_data$qvalue
          ylab_title <- "1 - qvalue"
        }else{
          if(method != "anova"){
            use_data %<>% .[.$P.adj < 0.05, ]
            use_data$Value <- 1 - use_data$P.adj
            ylab_title <- "1 - P.adjust"
          }else{
            stop("This function can not be used to ", method," currently!")
          }
        }
      }
    }
    if(!method %in% c("lefse", "rf")){
      if(length(unique(use_data$Comparison)) > 1){
        # make sure the Group not replicated for multiple comparisions
        if(is.null(select_group)){
          message('Multiple comparisions found. But select_group parameter not provided. Select the first group pair to show ...')
          select_group <- unique(use_data$Comparison)[1]
        }else{
          if(length(select_group) > 1){
            stop("The select_group parameter should only have one element! Please check the input!")
          }
          if(! select_group %in% use_data$Comparison){
            stop("The select_group parameter must be one of elements of object$res_diff$Comparison!")
          }
        }
        use_data %<>% .[.$Comparison %in% select_group, ]
      }
    }
    use_data %<>% .[!duplicated(.$Taxa), ]
    if(nrow(use_data) == 0){
      stop("No available data can be used to show!")
    }
    if(is.null(threshold)){
      sel_num <- use_number
    }else{
      sel_num <- sum(use_data$Value > threshold)
      if(sel_num == 0){
        stop("Too large threshold provided, no data selected!")
      }
      sel_num <- 1:sel_num
    }
    if(length(sel_num) > nrow(use_data)){
      sel_num <- 1:nrow(use_data)
    }
    use_data %<>% .[sel_num, ]
    if("group" %in% colnames(use_data)){
      if(is.null(group_order)){
        if((!is.null(datat$group_order)) & (length(unique(use_data$group)) == length(datat$group_order))){
          use_data$group %<>% factor(., levels = datat$group_order)
        }else{
          use_data$group %<>% as.character %>% as.factor
        }
      }else{
        use_data$group %<>% factor(., levels = group_order)
      }
      if(color_group_map){
        # fix colors for each group
        all_groups <- datat$group_order
        # make sure colors length enough for selection
        color_values <- expand_colors(color_values, length(all_groups))
        use_groups <- levels(use_data$group)
        color_values %<>% .[match(use_groups, all_groups)]
      }else{
        color_values <- expand_colors(color_values, length(levels(use_data$group)))
      }
      # rearrange orders
      if(length(levels(use_data$group)) == 2){
        use_data$Taxa %<>% as.character %>% factor(., levels = rev(unique(unlist(lapply(levels(use_data$group), function(x){
          if(x == levels(use_data$group)[1]){
            use_data[as.character(use_data$group) %in% x, ] %>% .[order(.$Value, decreasing = TRUE), "Taxa"]
          }else{
            use_data[as.character(use_data$group) %in% x, ] %>% .[order(.$Value, decreasing = FALSE), "Taxa"]
          }
        })))))
        use_data[use_data$group == levels(use_data$group)[2], "Value"] %<>% {. * -1}
      }else{
        use_data$Taxa %<>% as.character %>% factor(., levels = rev(unique(unlist(lapply(levels(use_data$group), function(x){
          use_data[as.character(use_data$group) %in% x, ] %>% .[order(.$Value, decreasing = TRUE), "Taxa"]
        })))))
      }
    }else{
      use_data %<>% .[order(.$Value, decreasing = TRUE), ]
      if(coord_flip){
        use_data$Taxa %<>% factor(., levels = rev(.))
      }else{
        use_data$Taxa %<>% factor(., levels = .)
      }
      ylab_title <- "Value"
    }
    datat$plot_diff_bar_taxa <- levels(use_data$Taxa) %>% rev

    if("group" %in% colnames(use_data)){
      p <- ggplot(use_data, aes(x = Taxa, y = Value, color = group, fill = group, group = group)) +
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values)
    }else{
      p <- ggplot(use_data, aes(x = Taxa, y = Value))
    }
    p <- p +
      geom_bar(stat = "identity", position = position_dodge(), ...) +
      theme_bw() +
      ylab(ylab_title) +
      xlab("") +
      theme(axis.title = element_text(size = 17), axis.text.y = element_text(size = axis_text_y, color = "black")) +
      theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank()) +
      theme(panel.border = element_blank()) +
      theme(axis.line.x = element_line(color = "grey60", linetype = "solid", lineend = "square"))

    if(coord_flip){
      p <- p + coord_flip()
    }else{
      p <- p + ggplot_xtext_anglesize(xtext_angle = xtext_angle, xtext_size = xtext_size)
    }
    p<-p+theme(text = element_text(family = "serif"),aspect.ratio = 0.5)
    ggsave(paste(export_path,"/microchat_lefse bar",".pdf",sep = ""),
           units = "cm",
           width = 21/3*2,
           height =  21*p$theme$aspect.ratio/3,
           p)
return(p)
}

tidy_dataset = function(self,main_data = FALSE){
  self <- tidy_samples(self)
  self$otu_table %<>% {.[apply(., 1, sum) > 0, , drop = FALSE]}
  taxa_list <- list(rownames(self$otu_table), rownames(self$tax_table), self$phylo_tree$tip.label) %>%
    .[!unlist(lapply(., is.null))]
  taxa_names <- Reduce(intersect, taxa_list)
  if(length(taxa_names) == 0){
    if(is.null(self$phylo_tree)){
      stop("No same feature names found between rownames of otu_table and rownames of tax_table! Please check rownames of those tables!")
    }else{
      stop("No same feature name found among otu_table, tax_table and phylo_tree! Please check feature names in those objects!")
    }
  }
  self$otu_table %<>% .[taxa_names, , drop = FALSE]
  if(!is.null(self$tax_table)){
    self$tax_table %<>% .[taxa_names, , drop = FALSE]
  }
  if(!is.null(self$phylo_tree)){
    self$phylo_tree %<>% ape::drop.tip(., base::setdiff(.$tip.label, taxa_names))
  }
  if(!is.null(self$rep_fasta)){
    if(!all(taxa_names %in% names(self$rep_fasta))){
      stop("Some feature names are not found in the names of rep_fasta! Please provide a complete fasta file or manually check the names!")
    }
    self$rep_fasta %<>% .[taxa_names]
  }
  # check again whether a sample has 0 abundance after feature filtering
  self <- tidy_samples(self)
  if(!main_data){
    sample_names <- rownames(self$sample_table)
    if(!is.null(self$taxa_abund)){
      self$taxa_abund %<>% lapply(., function(x) x[, sample_names, drop = FALSE])
    }
    if(!is.null(self$alpha_diversity)){
      self$alpha_diversity %<>% .[sample_names, , drop = FALSE]
    }
    if(!is.null(self$beta_diversity)){
      self$beta_diversity %<>% lapply(., function(x) x[sample_names, sample_names, drop = FALSE])
    }
  }
}

tidy_samples = function(microtable_obj){
  microtable_obj$otu_table <- check_abund_table(microtable_obj$otu_table)
  sample_names <- intersect(rownames(microtable_obj$sample_table), colnames(microtable_obj$otu_table))
  if(length(sample_names) == 0){
    stop("No same sample name found between rownames of sample_table and colnames of otu_table! ",
         "Please first check whether the rownames of sample_table are sample names! Then check through the sample names of each table!")
  }
  # keep the sample order same with original sample table
  sample_names <- rownames(microtable_obj$sample_table) %>% .[. %in% sample_names]
  microtable_obj$sample_table %<>% .[sample_names, , drop = FALSE]
  microtable_obj$otu_table %<>% .[ , sample_names, drop = FALSE]
  microtable_obj
}

check_abund_table = function(otu_table){
  if(!all(sapply(otu_table, is.numeric))){
    stop("Some columns in otu_table are not numeric class! Please check the input data!")
  }
  if(any(apply(otu_table, 1, sum) == 0)){
    remove_num <- sum(apply(otu_table, 1, sum) == 0)
    message(remove_num, " taxa with 0 abundance are removed from the otu_table ...")
    otu_table %<>% .[apply(., 1, sum) > 0, , drop = FALSE]
  }
  if(any(apply(otu_table, 2, sum) == 0)){
    remove_num <- sum(apply(otu_table, 2, sum) == 0)
    message(remove_num, " samples with 0 abundance are removed from the otu_table ...")
    otu_table %<>% .[, apply(., 2, sum) > 0, drop = FALSE]
  }
  if(ncol(otu_table) == 0){
    stop("No available sample! Please check the data!")
  }
  if(nrow(otu_table) == 0){
    stop("No available taxon! Please check the data!")
  }
  otu_table
}

summarySE_inter<-function (usedata = NULL, measurevar, groupvars = NULL, na.rm = TRUE)
{
  length2 <- function(x, na.rm = TRUE) ifelse(na.rm, sum(!is.na(x)),
                                              length(x))
  datac <- usedata %>% dplyr::grouped_df(groupvars) %>% dplyr::summarise(N = length2(!!sym(measurevar),
                                                                                     na.rm = na.rm), Mean = mean(!!sym(measurevar), na.rm = na.rm),
                                                                         SD = stats::sd(!!sym(measurevar), na.rm = na.rm)) %>%
    as.data.frame
  datac$SE <- datac$SD/sqrt(datac$N)
  datac
}



###plotALLDiffTree=====
"plotALLDiffInGlobalTree" <- function(submchat,diff_otus,layout="fan",open.angle=0,
                              branch.length=c("branch.length","none"),
                              colors_group=colorCustom(4,pal = "gygn"),
                              color_taxa=colorCustom(40,pal = "gygn"),
                              export_path) {
  branch.length<-match.arg(branch.length)
  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  abun <- submchat$otu_table
  taxon <- submchat$taxon_table
  tree<-submchat$tree

  sel.group<-unique(substr(colnames(abun),start = 1,stop = 2))
  sel.order<-match(sel.group,all.group)
  colors_group<-colors_group[sel.order]

  tree_filter<-tree_filter_second(tree,taxon=taxon,abun=abun)

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

  ttreex<-fortify(tree)
  ttreex$groupx<-ttreex$group
  ttreex$groupx<-as.character(ttreex$groupx)
  ttreex$groupx[!is.na(match(ttreex$label,setdiff(unique(ttreex$label),diff_otus)))]<-"grey"
  ttreex$groupx<-factor(ttreex$groupx,levels = unique(as.character(ttreex$groupx)))
  color_groupx<-c(color_group,"grey")
  names(color_groupx)<-c(names(color_group),"grey")

  p<-ggtree(ttreex,
            size=ifelse(ttreex$groupx=="grey",0.2,0.8),
            layout = layout,open.angle = open.angle,
            aes(color=groupx),
            ladderize=TRUE,branch.length = branch.length)+
    scale_color_manual(values=color_groupx,
                       guide=guide_legend(keywidth = 1, name = "\nPhylum",
                                          keyheight = 1, order=1,
                                          override.aes=list(size=5,shape = 20)))+
    guides(color="none",size="none")

  p<-p+ new_scale_fill()+new_scale_color()+
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

  abun_sumx<-abun_sum1[which(abun_sum1$id %in% diff_otus),]
  p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sumx, geom=geom_point, ###geom_star
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
                                         override.aes=list(size=2)))+
    guides(color="none",fill="none")

  p <- p +
    #geom_treescale(x=5,y=5,fontsize=2, width=5,linesize=0.3,family="serif",color="grey50")+
    theme(legend.position=c(0, 0.8),
          text=element_text(family = "serif"),
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=7),
          legend.text=element_text(size=6),
          legend.spacing.y = unit(0.02, "cm")
    )

  export_pathx<-paste(export_path,"/data_microbiome/microbial differential analysis",sep = "")
  ggsave(paste(export_pathx,"/","AllDiffOTU_in_GlobalTree.pdf",sep = ""),p,
         width = 21,height = 21/2,units = "cm")
  #cat("please seek ",paste(export_pathx,"/","AllDiffOTU_in_GlobalTree.pdf",sep = ""))
  return(p)
}

"plotALLDiffTree" <- function(submchat,diff_otus,specificOTU,layout="fan",open.angle=0,
                              branch.length=c("branch.length","none"),
                              colors_group=colorCustom(4,pal = "gygn"),
                              color_taxa=colorCustom(40,pal = "gygn"),
                              export_path) {
  branch.length<-match.arg(branch.length)
  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  abun <- submchat$otu_table
  taxon <- submchat$taxon_table
  tree<-submchat$tree

  sel.group<-unique(substr(colnames(abun),start = 1,stop = 2))
  sel.order<-match(sel.group,all.group)
  colors_group<-colors_group[sel.order]

  tree_2<-treeio::drop.tip(tree,setdiff(tree$tip.label,diff_otus))
  tree_filter<-tree_filter_second(tree_2,taxon=taxon,abun=abun)

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


  ttreex<-fortify(ttree)
  ttreex$groupx<-ttreex$group
  ttreex$groupx<-as.character(ttreex$groupx)
  ttreex$groupx[!is.na(match(ttreex$label,setdiff(unique(ttreex$label),specificOTU)))]<-"grey"
  ttreex$groupx<-factor(ttreex$groupx,levels = unique(as.character(ttreex$groupx)))
  color_groupx<-c(color_group,"grey")
  names(color_groupx)<-c(names(color_group),"grey")


  p<-ggtree(ttreex,
            size=ifelse(ttreex$groupx=="grey",0.3,0.8),
            layout = layout,open.angle = open.angle,
            aes(color=group),
            ladderize=TRUE,branch.length = branch.length)+
    scale_color_manual(values=color_group,
                       guide=guide_legend(keywidth = 1, name = "\nPhylum",
                                          keyheight = 1, order=1,
                                          override.aes=list(size=5,shape = 20)))+
    guides(color="none",size="none")

  abun_sum2<-abun_sum1[which(abun_sum1$id%in%diff_otus),]
  p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum2, geom=geom_bar, ###geom_star
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


  abun_sum2x<-abun_sum2[which(abun_sum2$id %in% diff_otus1),]
  p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum2x, geom=geom_point, ###geom_star
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
                                         override.aes=list(size=2)))+
    guides(color="none",fill="none")

  p <- p +
    #geom_treescale(x=5,y=5,fontsize=2, width=5,linesize=0.3,family="serif",color="grey50")+
    theme(legend.position=c(0, 0.8),
          text=element_text(family = "serif"),
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=7),
          legend.text=element_text(size=6),
          legend.spacing.y = unit(0.02, "cm")
    )

  export_pathx<-paste(export_path,"/data_microbiome/microbial differential analysis",sep = "")
  ggsave(paste(export_pathx,"/","AllDiffOTU_Tree.pdf",sep = ""),p,
         width = 21,height = 21/2,units = "cm")
  #cat("please seek ",paste(export_pathx,"/","AllDiffOTU_Tree.pdf",sep = ""))
  return(list(plot=p,data=ttreex,color=color_groupx))
}

'plotSpecificDiffTree' <- function(aLlDifftree,specificOTU,export_path) {
  ttreex<-aLlDifftree$data
  ttreexx<-ttreex[which(ttreex$label %in%specificOTU),]
  p<-aLlDifftree$plot
  color_groupx<-aLlDifftree$color

  p+geom_tiplab(data=ttreexx, hjust=-1,
                mapping=aes(color=group),
                align = TRUE, linetype = "dotted",
                linesize = 0.5,size=0)+
    scale_color_manual(values = color_groupx)->p

  export_pathx<-paste(export_path,"/data_microbiome/microbial differential analysis",sep = "")
  ggsave(paste(export_pathx,"/","SpecificDiffOTU_Tree.pdf",sep = ""),p,
         width = 21,height = 21/2,units = "cm")
  #cat("please seek ",paste(export_pathx,"/","SpecificDiffOTU_Tree.pdf",sep = ""))
  return(p)
}

