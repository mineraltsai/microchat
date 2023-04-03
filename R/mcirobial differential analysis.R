"calcMicrochatLefse" <- function(submchat,
                               lda_score=2,
                               file.save=TRUE,
                               export_path="microbial differential analysis") {
  dir.create(export_path, recursive = TRUE)
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

  file2=paste(export_path,"/lefse result (lda: ",lda_score,") .txt",sep = "")
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
 suppressMessages(library(ggplot2))
  dir.create(export_path, recursive = TRUE)
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
                        alpha=0.5,
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
                          alpha=0.5,
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

  if (file.save) ggsave(paste(export_path,"/lefse cladogram (with legend) of microbial differential analysis",".pdf",sep = ""),p_legend)
  if (file.save) ggsave(paste(export_path,"/lefse cladogram (without legend)of microbial differential analysis",".pdf",sep = ""),p_nolegend)

  cat("\n","Lefse cladogram of differential analysis has been exported. Please check it.","\n")

  return(list(p_nolegend=p_nolegend,p_legend=p_legend))
}


"plotMicrochatComplexLefse" <- function(submchat,lda_score=2,
                                        control="ct",
                                        layout="radial",
                                        label_plot_display=4,
                                        color_group=c("#FF0000","#333399","#009933","#00FFFF","#EC00FF","yellow"),
                                        export_path="microbial differential analysis") {

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
                          color_group=color_group,
                          label_plot_display=label_plot_display,
                          file.save=FALSE,
                          export_path=export_path)

    ggsave(paste(export_path,"/lefse_nolegend (",control,"-",kk,")",".pdf",sep = ""),p$p_nolegend)
    ggsave(paste(export_path,"/lefse_legend (",control,"-",kk,")",".pdf",sep = ""),p$p_legend)
    cat("Pairwise lefse cladogram have been exported to","/",export_path,"",sep = "")
  }
}



"calcMicrochatLimma" <- function(submchat,
                                 comparison = "ml-ct",
                                 file.save=TRUE,
                                 export_path="microbial differential analysis/limma") {
  message(" message: ajust.method=BH,test=glm")
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

  select_comparison<-strsplit(comparison,"-")
  select_comparison<-select_comparison[[1]]

  select_comparison<-paste(select_comparison[2],
                           select_comparison[1],sep = "-")

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

  file2=paste(export_path,"/Differential analysis results (limma: ",select_comparison,") .txt",sep = "")
  if (file.save) write.table(result,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Differential analysis results based on edgeR and limma (",select_comparison,") has been exported to ","/",export_path,"",sep = "","\n")

  microchatLimmaobj<-list(result,
                          select_comparison)

  names(microchatLimmaobj)[1]<-"result"
  names(microchatLimmaobj)[2]<-"select_comparison"

  class(microchatLimmaobj) <- c("microchat","data.frame")

  return(microchatLimmaobj)
}






"plotMicrochatDiffVolcano" <- function(microchatLimmaobj,
                                       sigcolor=c('#a121f0','#bebebe','#ffad21'),
                                       sigalpha=c(1,0.4,1),
                                       layout="choice1",
                                       annosize=3,
                                       pvalue_thres=0.05,
                                       foldchange_thres=1,
                                       export_path="difference analysis/limma"

) {
  dir.create(export_path, recursive = TRUE)

  if (class(microchatLimmaobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data<-microchatLimmaobj
  select_comparison<-data$select_comparison
  data<-data$result
  data %>%
    mutate(change = case_when(
      pvalue < pvalue_thres & logFC > foldchange_thres ~ "UP",
      pvalue < pvalue_thres & logFC < -foldchange_thres ~ "DOWN",
      TRUE ~ "NOT"
    )) -> DEG

  label<-table(DEG$change)


  if (layout=="choice1") {
    ggplot(data=DEG,aes(x=logFC,
                        y=-log10(pvalue),

                        color=change))+
      geom_point(alpha=0.8,size=3)+
      labs(x="log2 fold change")+ ylab("-log10 pvalue")+
      #ggtitle(this_title)+
      theme_bw(base_size = 20)+
      #theme(plot.title = element_text(size=15,hjust=0.5),)+
      scale_color_manual(values=sigcolor) -> p1

    p1 +
      geom_vline(xintercept = foldchange_thres,lty="dashed")+
      geom_vline(xintercept = -foldchange_thres,lty="dashed") -> p2
    p2<-p2+
      annotate('text', size = annosize,color=sigcolor[1],
               family="serif",label = paste("Down: ",label[1],sep = ""), min(DEG$logFC)*0.9, max((-log10(data$pvalue)))*0.9) +
      annotate('text',size = annosize,color=sigcolor[3],
               family="serif", label = paste("Up: ",label[3],sep = ""), max(DEG$logFC)*0.9, max((-log10(data$pvalue)))*0.9) +
      theme(aspect.ratio = 1,
            axis.title.x = element_text(family = "serif"),
            axis.title.y = element_text(family = "serif")
      )

  }


  if (layout=="choice2"){
    ####choice2
    table(DEG$change)
    otu2<-DEG[which(DEG$change != "NOT"),]
    otu2$name<-rownames(otu2)
    otu2<-otu2[1:10,]
    otu<-DEG
    otu[which(otu$change=="UP" |otu$change=="DOWN"),"sigalpha"]<-"sig"
    otu[which(otu$change!="UP" & otu$change!="DOWN"),"sigalpha"]<-"notsig"
    table(otu$sigalpha)

    p2<-ggplot(otu, aes(x = logCPM, y = logFC,
                        color = change,
                        alpha= change,
                        size=-log10(otu$pvalue))) +
      geom_point() +
      scale_alpha_manual(values = sigalpha,guide="none") +
      scale_color_manual(values = sigcolor,
                         guide = "none") +
      theme(legend.position = "none")+
      ggrepel::geom_text_repel(data=otu2,
                               aes(x = logCPM, y = logFC, label = name),
                               size = 2.5,
                               box.padding = unit(0.3, 'lines'), show.legend = FALSE,
                               max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
      theme(aspect.ratio = 1) +
      annotate('text', size = annosize,color=sigcolor[3],
               family="serif",label = paste("Up: ",label[3],sep = ""), max(DEG$logCPM)*0.9, max(((data$logFC)))*0.9) +
      annotate('text',size = annosize,color=sigcolor[1],
               family="serif", label = paste("Down: ",label[1],sep = ""), max(DEG$logCPM)*0.9, min(((data$logFC)))*0.9) +
      theme(panel.background=element_rect(fill="grey95"),
            legend.position="right",
            plot.margin=margin(0,0,0,0),
            axis.title.x=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = -1),
            axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
            legend.spacing.y=unit(0.02, "cm"),
            legend.title=element_text(size=12,face = "bold",family = "serif"),
            legend.text=element_text(size=8,face = "bold",family = "serif"),
            legend.box.spacing=unit(0.02,"cm"))


  }
  p2

  ggsave(paste(export_path,"/single_volcano_",paste(select_comparison[2],select_comparison[1],sep = "-"),".pdf",sep = ""),p2)
  cat("Single volcano has been exported to","/",export_path,"",sep = "")

  file2=paste(export_path, "/single_volcano_",paste(select_comparison[2],select_comparison[1],sep = "-"),".txt",sep = "" )
  write.table(data,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  cat("\n","Differential analysis results has been exported to“",export_path,"",sep = "")

  return(list(data=DEG,plot=p2,sig_condtion=label))
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


  file2=paste(export_path,"/Muti-differential analysis results (limma: ",comparison,") .txt",sep = "")
  if (file.save) write.table(cpma,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Muti-differential analysis results based on edgeR and limma (",select_comparison,") has been exported to ","/",export_path,"",sep = "","\n")

  microchatComplexLimmaobj<-list(data=cpma,
                                 select_comparison1=select_comparison1,
                                 select_comparison2=select_comparison2)


  names(microchatComplexLimmaobj)[1]<-"data"
  names(microchatComplexLimmaobj)[2]<-"select_comparison1"
  names(microchatComplexLimmaobj)[3]<-"select_comparison2"

  class(microchatComplexLimmaobj) <- c("microchat","data.frame")

  return(microchatComplexLimmaobj)
}


"plotMicrochatComplexDiffVolcano" <- function(microchatComplexLimmaobj,
                                              repelsize=3,
                                              repeldisp_rate=2, #1-all,2-half
                                              annosize=3,
                                              pvalue=0.05,
                                              foldchange=1,
                                              sigcolor=c('red', 'orange', 'green', 'blue'),
                                              export_path="difference analysis/limma") {
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
          axis.title.x=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = -1),
          axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
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
  p2 <- p1 + theme(legend.position = 'right') +
    ggrepel::geom_text_repel(data = rbind(up, down),
                             aes(x = logFC_1, y = logFC_2, label = rownames),
                             size = repelsize,box.padding = unit(0.5, 'lines'), family="serif",
                             segment.color = 'black', show.legend = FALSE)


  ggsave(paste(export_path,"/Muti-volcano (",control,"-",treat1,"-",treat2,").pdf",sep = ""),p2)
  cat("Muti-volcano plot has been exported to","/",export_path,"",sep = "")
  return(p2)

}
