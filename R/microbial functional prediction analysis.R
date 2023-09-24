"plot_tax4fun_extendbar" <- function(submchat,
                                     padj=FALSE,
                                     comparison_sel="ct-xl",
                                     ko1.select="Metabolism",   ###  ko1.select="all",则使用所有ko1
                                     min.fun.abun = 0.3,   ###平均丰度大于该值才考虑组图
                                     filter=FALSE,          ###是否过滤掉差异不显著的
                                     layout.rt=c(6,4,4,1.5),
                                     add.bar.color=colorCustom(50,pal = "gygn"),
                                     method="t.test",  ### "wilcox.test"
                                     color_group=c("#FF0000","#333399","#009933","#00FFFF","#EC00FF","yellow"),
                                     silva_path='/Volumes/BOOTCAMP/SILVA123',
                                     export_path="tax4fun") {

  export_path<-paste(export_path,"/microbial functional prediction analysis",sep = "")
  dir.create(export_path, recursive = TRUE)

  library(Tax4Fun)
  data(kegg_anno)
  kegg_anno=kegg_anno

  #all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))

  abun<-submchat$otu_table
  #sel.group<-unique(substr(colnames(abun),start = 1,stop = 2))
  tt<-str_split_fixed(comparison_sel,'-',2)%>%as.character()

  #sel.order<-match(sel.group,all.group)
  #color_group<-color_group[sel.order]

  color.all<-color_group
  allg<-unique(substr(colnames(abun),start = 1,stop = 2))
  color.use<-color.all[match(tt,allg)]
  #color.use<-color_group

  qiimeotu<-  readtaxdata(submchat,comparison = comparison_sel)
  QIIMESingleData<-tax4funpretreat(qiimeotu)
  tax4fp<-tax4funpred(QIIMESingleData,
                      kegg_anno,
                      folderReferenceData = silva_path,
                      export_path=export_path)
  group<-qiimeotu$sampleda
  fungene<-tax4fp$tax4fun_gene%>%rownames_to_column(var = "fungene")
  funpred<-tax4fp$tax4fun_pathway

  library(cowplot)
  pp<-plot_funpredict(funpred,
                      padj=padj,
                      sampledata = qiimeotu$sampleda,
                      ko1.select=ko1.select,   ###  ko1.select="all",则使用所有ko1
                      min.fun.abun = min.fun.abun,   ###平均丰度大于该值才考虑组图
                      filter=filter,          ###是否过滤掉差异不显著的
                      sig.text.color="black",
                      sig.text.size=2,
                      ko2.text.color="black",
                      ko2.text.size=4,
                      add.bar.text.color="white",
                      add.bar.text.size=3,
                      layout.rt=layout.rt,
                      add.bar.color=add.bar.color,
                      compare.color=color.use,
                      method=method,  ### "wilcox.test"
                      siglabel="sig+label") ### "sig"  "label",
  pp
  id<-group$group%>%unique()
  ggsave(paste(export_path,"/",paste(id[1],id[2],sep = "-"),".pdf",sep = ""),
         pp)
  cat("图片已保存至相对路径","/",export_path,"","下",sep = "")
  return(list(extendbar=pp,pred_res=funpred))
}



"readtaxdata" <- function(submchat,comparison = "M0-M5")
{
  comparison_sel<-comparison
  data_prep_new<-data_prep(submchat)
  data<-data_prep_new$data
  group_order<-data_prep_new$group


  taxdata<-submchat$taxon_table

  taxon1<-taxdata
  for (k in 1:ncol(taxdata)) {
    taxon1[,k]<-substr(taxdata[,k],start = 4, stop = 100)%>%data.frame()
    colnames(taxon1)[k]<-colnames(taxdata)[k]
    rownames(taxon1)<-rownames(taxdata)
  }
  taxon<-taxon1

  counts<-data$otu_table
  sampledata<-data$sampledata
  sampleda<-rownames_to_column(sampledata,var = "rowname")
  colnames(sampleda)[1]<-"sample"

  if (comparison_sel!="all") {
    select_comparison<-strsplit(comparison_sel,"-")
    select_comparison<-select_comparison[[1]]
  } else {
    select_comparison<-sampleda$group
  }

  slee<-match(sampleda$group,select_comparison)%>%as.numeric()
  pos<-which(!is.na(slee))

  abun_sel<-submchat$otu_table[,pos]
  sampleda<-sampleda[pos,]

  taxonomy<-unite(taxon,"taxonomy",sep = "; ")
  taxonomy<-rownames_to_column(taxonomy,var = "geneid")
  abun_sel<-rownames_to_column(abun_sel,var = "geneid")

  abun_sel1<-merge(abun_sel,taxonomy,by="geneid")

  col_names <- colnames(abun_sel1)
  col_classes <- rep("numeric", times = length(col_names))
  col_classes[1] <- "character"

  col_classes[length(col_classes)] <- "character"

  full_otu_table <- abun_sel1

  data_cols <- 2:(length(col_names) - 1)

  sample_ids <- col_names[data_cols]
  otu_ids <- as.character(full_otu_table[, 1])
  counts <- as.matrix(full_otu_table[, data_cols])
  rownames(counts) <- otu_ids

  metadata_vals <- as.character(full_otu_table[, length(col_names)])
  names(metadata_vals) <- otu_ids

  tt<-list(sample_ids = sample_ids, otu_ids = otu_ids, counts = counts,
           sampleda=sampleda,treat=select_comparison[1],
           metadata = metadata_vals)

  return(tt)

}


"tax4funpretreat" <- function(qiimeotu) {

  ModSilvaIds <- gsub("uncultured archaeon", "", qiimeotu$metadata)
  ModSilvaIds <- gsub("uncultured organism", "", ModSilvaIds)
  ModSilvaIds <- gsub("uncultured bacterium", "", ModSilvaIds)
  ModSilvaIds <- gsub("uncultured crenarchaeote", "", ModSilvaIds)
  ModSilvaIds <- gsub("uncultured euryarchaeote", "", ModSilvaIds)
  ModSilvaIds <- gsub("; ", ";", ModSilvaIds)
  otuTable <- rowsum(data.frame(qiimeotu$counts), ModSilvaIds)
  inputData <- list(sampleNames = qiimeotu$sample_ids,
                    otuTable = otuTable)

  return(inputData)
}

"tax4funpred" <- function(QIIMESingleData=QIIMESingleData,
                          kegg_anno=kegg_anno,
                          folderReferenceData = '/Volumes/BOOTCAMP/SILVA123',
                          export_path='~/Desktop/文章撰写/jbh1/tax4fun') {

  #KEGG 功能基因（KO 第四级）丰度预测
  Tax4FunOutput <- tax4Fun(QIIMESingleData,
                           folderReferenceData,
                           fctProfiling = TRUE,
                           refProfile = 'UProC',
                           shortReadMode = TRUE,
                           normCopyNo = TRUE)
  tax4fun_gene <- as.data.frame(t(Tax4FunOutput$Tax4FunProfile))
  gene <- rownames(tax4fun_gene)
  groupname<-unique(substr(colnames(tax4fun_gene),start = 1,stop = 2))
  groupname.use<-paste(groupname[1],groupname[2],sep = "_")

  cat("\n",paste("文件位于'",export_path,"/' 下",sep = ""),"\n")
  dir.create(export_path, recursive = TRUE)  #创建目录“plot”用于存放图片
  write.table(cbind(gene, tax4fun_gene),
              paste(export_path,"/fungene.",groupname.use,".txt",sep = ""),
              row.names = FALSE, sep = '\t', quote = FALSE)

  #KEGG 代谢通路（KO 第三级）丰度预测
  Tax4FunOutput <-tax4Fun(QIIMESingleData,
                          folderReferenceData,
                          fctProfiling = FALSE,
                          refProfile = 'UProC',
                          shortReadMode = TRUE,
                          normCopyNo = TRUE)
  tax4fun_pathway <- as.data.frame(t(Tax4FunOutput$Tax4FunProfile))
  pathway <- rownames(tax4fun_pathway)
  write.table(cbind(pathway, tax4fun_pathway),
              paste(export_path,"/pathway.",groupname.use,".txt",sep = ""),
              row.names = FALSE, sep = '\t', quote = FALSE)

  ##为 KEGG 代谢通路映射 KO 分类
  tax4fun_pathway$KO3_id <- t(data.frame(strsplit(rownames(tax4fun_pathway), ';')))[, 1]
  tax4fun_pathway$KO3_id <- gsub('ko', '', tax4fun_pathway$KO3_id)
  tax4fun_pathway <- merge(tax4fun_pathway, kegg_anno, by = 'KO3_id')
  write.table(tax4fun_pathway,
              paste(export_path,"/pathway.anno.",groupname.use,".txt",sep = ""),
              row.names = FALSE, sep = '\t', quote = FALSE)

  return(list(tax4fun_gene=tax4fun_gene,tax4fun_pathway=tax4fun_pathway))
}


"plot_funpredict" <- function(funpred,
                              sampledata,
                              padj=FALSE,
                              ko1.select="Metabolism",
                              min.fun.abun = 0.3,
                              filter=TRUE,
                              sig.text.size=3,
                              sig.text.color="red",
                              ko2.text.size=2,
                              ko2.text.color="red",
                              add.bar.text.size=3,
                              add.bar.text.color="white",
                              layout.rt=c(4,4,4,1.5),
                              add.bar.color=colorCustom(50,pal = "gygn"),
                              compare.color=colorCustom(2,pal = "gygn"),
                              method=c("wilcox.test","t.test"),
                              siglabel=c("sig","label","sig+label")) {

  method<-match.arg(method)
  siglabel<-match.arg(siglabel)

  group<-sampledata
  a<-subset(funpred,select=c(2:(nrow(group)+1),KO1,KO2,KO3))
  datako2<-subset(funpred,select=c(2:(nrow(group)+1),KO2))
  datako3<-subset(funpred,select=c(2:(nrow(group)+1),KO3))
  datako4<-subset(funpred,select=c(2:(nrow(group)+1),KO1,KO3))
  datako4<-unite(datako4,"KO3",KO1,KO3,sep = ";")

  datako4 %>%
    group_by(KO3) %>%
    summarise_all(sum) -> data
  data[,c("KO1","KO3")]<-str_split_fixed(data$KO3,pattern=";",2)

  if (ko1.select !="all") {
    data<-data[which(data$KO1==ko1.select),]
  } else {
    data<-data
  }


  data<-subset(data,select = c(-KO1))
  data<-column_to_rownames(data,var = "KO3")

  data <- data*100
  data <- data %>% filter(apply(data,1,mean) > min.fun.abun)
  data <- t(data)
  data1 <- data.frame(data,group$group)
  colnames(data1) <- c(colnames(data),"Group")
  data1$Group <- as.factor(data1$Group)

  if (method=="t.test") {  ## t-test

    diff <- data1 %>%
      select_if(is.numeric) %>%
      map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')

    if (padj) diff$p.value <- p.adjust(diff$p.value,"fdr")
    if (filter) {
      diff1 <- diff %>% filter(p.value < 0.05)
      if (nrow(diff1) <= 0) {
        diff<-diff
      } else {
        diff<-diff1
      }
    }
  }

  if (method=="wilcox.test"){  ## wilcox
    library(tidyverse)
    diff <- data1 %>%
      select_if(is.numeric) %>%
      map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1)), .id = 'var')

    if (padj) diff$p.value <- p.adjust(diff$p.value,"fdr")
    if (filter) {
      diff1 <- diff %>% filter(p.value < 0.05)
      if (nrow(diff1) <= 0) {
        diff<-diff
      } else {
        diff<-diff1
      }
    }
  }


  ## 绘图数据构建
  ## 左侧条形图
  abun.bar <- data1[,c(diff$var,"Group")] %>%
    gather(variable,value,-Group) %>%
    group_by(variable,Group) %>%
    summarise(Mean = mean(value))

  kk1<-subset(funpred,select=c(KO1,KO2,KO3))
  kk1$variable<-kk1$KO3
  KKK<-merge(abun.bar,kk1,by="variable")
  abun.bar<-KKK



  ## 右侧散点图
  diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
  diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                              levels(data1$Group)[2]))
  diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]

  ## 左侧条形图
  library(ggplot2)
  cbbPalette <- compare.color


  abun.bar<-abun.bar[order(abun.bar$KO2),]
  rownames(abun.bar)<-1:nrow(abun.bar)
  abun.bar$variable<-factor(abun.bar$variable,levels=rev(unique(abun.bar$variable)))
  names(cbbPalette)<-unique(group$group)


  p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
    scale_x_discrete(limits = levels(abun.bar$variable),label=NULL) +
    coord_flip() +
    xlab("") +
    ylab("Mean proportion (%)") +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length.y = unit(0,"lines"),
          axis.ticks.y  = element_line(size=0),
          axis.line  = element_line(colour = "black"),
          axis.title.y=element_text(colour='black', size=8,face = "bold"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          legend.title=element_blank(),
          legend.background = element_blank(),
          legend.text=element_text(size=8,face = "bold",colour = "black",
                                   family = "serif",
                                   margin = margin(r = 20)),
          legend.position = c(-0.5,-0.05),
          legend.direction = "horizontal",
          legend.key = element_rect(colour = "white"),
          legend.key.width = unit(0.8,"cm"),
          legend.key.height = unit(0.5,"cm"))

  for (i in 1:(nrow(diff.mean) - 1)) {
    p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }

  diff.mean$KO3<-diff.mean$var
  ajj<-merge(diff.mean,kk1,by="KO3", all = TRUE,sort = FALSE)
  df<-filter(ajj,ajj$var!="NA")

  color <- c( '#984EA3', '#FFED6F','#FFFF33','#377EB8',
              '#00D098', '#CCEBC5', '#4DAF4A','#E41A1C','#FDB462',
              '#B3DE69', '#BC80BD', '#40ffd1','#B0FF39','#FCCDE5',
              '#B02664', '#31a2ff','#ff0000','#76c8ff')
  names(color) <- names(df$KO2)

  p1 <- p1 +
    geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
    scale_fill_manual(values=cbbPalette)


  colorsbar<-add.bar.color[1:length(unique(abun.bar$KO2))]
  names(colorsbar)<-unique(abun.bar$KO2)
  abun.bar$KO4<-abun.bar$KO2
  abun.bar$KO4[setdiff(1:nrow(abun.bar),match(unique(abun.bar$KO4),abun.bar$KO4))]<-""
  abun.bar$KO4[match(unique(abun.bar$KO4),abun.bar$KO4)]<-abun.bar$KO4[match(unique(abun.bar$KO4),abun.bar$KO4)]
  # need to set `coord_flip = TRUE` if you plan to use `coord_flip()`
  ydens <- axis_canvas(p1, axis = "y", coord_flip = TRUE) +
    scale_x_discrete()+
    coord_flip()+
    geom_col(data=abun.bar, aes(x=variable,
                                y=1,
                                fill=KO2), alpha=1, width = 1,
             size=.2) +
    geom_text(data=abun.bar, aes(x=variable,
                                 y=2,
                                 label=variable),size=add.bar.text.size,
              color=add.bar.text.color,
              family="serif",check_overlap = FALSE,
              hjust = 1)+
    geom_text(data=abun.bar, aes(x=variable,
                                 y=0,
                                 label=KO4),size=ko2.text.size,
              color=ko2.text.color,
              family="serif",check_overlap = FALSE,
              hjust = 0)+
    scale_fill_manual(values = colorsbar)



  ## 右侧散点图
  abun.bar33<-subset(abun.bar,select=c(KO2,KO3))
  abun.bar33<-abun.bar33[duplicated(abun.bar33$KO3),]
  rownames(abun.bar33)<-1:nrow(abun.bar33)

  diff.mean[,"KO2"]<-abun.bar33$KO2[match(diff.mean$KO3,abun.bar33$KO3)]
  diff.mean<-diff.mean[order(diff.mean$KO2,diff.mean$KO3),]

  diff.mean$KO2 <- factor(diff.mean$KO2,levels = rev(unique(diff.mean$KO2)))
  diff.mean$p.value <- signif(diff.mean$p.value,3)
  diff.mean$p.value <- as.character(diff.mean$p.value)
  diff.mean$var<-factor(diff.mean$var,levels =rev(diff.mean$var))

  diff.mean[,"sig"]<-ifelse(diff.mean$p.value<0.001,"***",
                            ifelse(diff.mean$p.value<0.01,"**",
                                   ifelse(diff.mean$p.value<0.05,"*","ns") ))


  p2 <- ggplot(diff.mean,aes(x=var,y=estimate,fill = Group)) +
    scale_x_discrete(limits = levels(abun.bar$variable)) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +

    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals")


  for (i in 1:(nrow(diff.mean) - 1))
    p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

  p2 <- p2 +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(0.8), width = 0.5, size = 0.5) +
    geom_point(shape = 21,size = 3) +
    scale_fill_manual(values=cbbPalette) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')



  if (siglabel=="label"){
    p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
      geom_text(aes(y = 0,x = var),label = diff.mean$p.value,family="serif",
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = sig.text.size) +
      geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
                family="serif",srt = 90,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_text(family = "serif"),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }

  if (siglabel=="sig"){
    p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
      geom_text(aes(y = 0,x = var),label = diff.mean$sig,family="serif",
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = sig.text.size) +
      geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
                family="serif",srt = 90,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_text(family = "serif"),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }

  if (siglabel=="sig+label"){
    p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
      scale_x_discrete(limits = levels(abun.bar$variable)) +
      geom_text(aes(y = 0,x = var),label = paste(diff.mean$p.value,diff.mean$sig,sep = ""),
                family="serif",color=sig.text.color,
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = sig.text.size) +
      geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
                family="serif",srt = 90,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_text(family = "serif"),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }

  ## 图像拼接
  library(patchwork)


  p<-ydens+p1 + p2 + p3 + plot_layout(widths = c(layout.rt[1],
                                                 layout.rt[2],
                                                 layout.rt[3],
                                                 layout.rt[4]),byrow=TRUE)

  cat("\n已输出",
      paste(unique(substr(unique(group$sample),start = 1,stop = 2))),"基于",method,"的Extended error bar plot (扩展柱状图)")
  return(p)

}



"plot_tax4fun_cladogram" <- function(submchat,
                                     ko1.select="Metabolism", ###如不为all,则可以为Metabolism，则仅可视化KO2-KO3两个层次;如为all，则可视化KO1-KO2-KO3三个层次
                                     ldascore=2,
                                     layout="radial",#布局类型
                                     color_group=c("#00AED7", "#FD9347","#ff0000","#ece12e","red"),
                                     silva_path='/Volumes/BOOTCAMP/SILVA123',
                                     export_path="tax4fun"
) {
  dir.create(export_path, recursive = TRUE)
  data(KO)
  KO = KO
  KO$name <- NULL
  KO <- as.data.frame(KO) %>%
    unnest(cols = c("children.name","children.children"),names_repair = tidyr_legacy) %>%#重要函数
    unnest(cols = c("children.name","name","children"),names_repair = tidyr_legacy) %>%
    unnest(cols = c("children.name","name","name1","children"),names_repair = tidyr_legacy)
  colnames(KO) <- c("L1","L2","L3","KO")
  KO %<>% #整理KEGG ORTHOLOGY
    separate(col = "KO",sep = "  ",into = c("KO","Description")) %>%
    separate(col = "L1",sep = " ",into = c("L1_ID","L1"),extra = "merge") %>%
    filter(!L1_ID %in% c("09180","09190")) %>% #去除BRITE hierarchies和Not Included in Pathway or Brite两大类
    separate(col = "L2",sep = " ",into = c("L2_ID","L2"),extra = "merge") %>%
    separate(col = "L3",sep = " ",into = c("L3_ID","L3"),extra = "merge") %>%
    separate(col = "L3",sep = " \\[PATH:",into = c("L3","PathwayID")) %>%
    mutate(PathwayID=str_remove(PathwayID,pattern = "\\]")) %>%
    drop_na()#KEGG ORTHOLOGY等级有缺失的删掉
  head(KO)

  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  abun<-submchat$otu_table
  sel.group<-unique(substr(colnames(abun),start = 1,stop = 2))

  sel.order<-match(sel.group,all.group)
  color_group<-color_group[sel.order]
  names(color_group)<-sel.group
  groupnum<-unique(substr(colnames(abun),start = 1,stop = 2))%>%length() ###组数量
  groupp<-unique(substr(colnames(abun),start = 1,stop = 2)) ###组名

  if (groupnum %% 2==0) {
    pair<-groupnum/2
    fungene_all<-data.frame(fungene=0)
    for (tt in 1:pair) {
      k<-tt*2-1
      sel_group<-groupp[k:(k+1)]
      sel_group1<-paste(sel_group[1],sel_group[2],sep = "-")

      qiimeotu<-  readtaxdata(submchat,comparison = sel_group1)

      QIIMESingleData<-tax4funpretreat(qiimeotu)

      data(kegg_anno)
      kegg_anno=kegg_anno
      tax4fp<-tax4funpred(QIIMESingleData,
                          kegg_anno,
                          folderReferenceData = silva_path,
                          export_path=export_path)

      fungene<-tax4fp$tax4fun_gene%>%rownames_to_column(var = "fungene")

      {
        fungene_all<-merge(fungene_all,fungene,by="fungene",all = TRUE)
      }
    }
    fungene_all<-fungene_all[-1,]
    rownames(fungene_all)<-1:nrow(fungene_all)
  } else {

    pair<-round(groupnum/2)
    fungene_all<-data.frame(fungene=0)
    for (tt in 1:(pair-1)) {
      k<-tt*2-1
      sel_group<-groupp[k:(k+1)]
      sel_group1<-paste(sel_group[1],sel_group[2],sep = "-")

      qiimeotu<-  readtaxdata(submchat,comparison = sel_group1)

      QIIMESingleData<-tax4funpretreat(qiimeotu)

      data(kegg_anno)
      kegg_anno=kegg_anno
      tax4fp<-tax4funpred(QIIMESingleData,
                          kegg_anno,
                          folderReferenceData = silva_path,
                          export_path=export_path)

      fungene<-tax4fp$tax4fun_gene%>%rownames_to_column(var = "fungene")

      {
        fungene_all<-merge(fungene_all,fungene,by="fungene",all = TRUE)
      }
    }
    fungene_all<-fungene_all[-1,]
    rownames(fungene_all)<-1:nrow(fungene_all)

    sel_group<-groupp[(pair*2-2):(pair*2-1)]
    sel_group1<-paste(sel_group[1],sel_group[2],sep = "-")

    qiimeotu<-  readtaxdata(submchat,comparison = sel_group1)

    QIIMESingleData<-tax4funpretreat(qiimeotu)

    data(kegg_anno)
    kegg_anno=kegg_anno

    tax4fp<-tax4funpred(QIIMESingleData,
                        kegg_anno,
                        folderReferenceData = silva_path,
                        export_path=export_path)

    fungene<-tax4fp$tax4fun_gene%>%rownames_to_column(var = "fungene")
    fungene_all<-merge(fungene_all,fungene,by="fungene",all = TRUE)
  }

  ffungene<-fungene_all
  ffungene[is.na(ffungene)]<-"0"
  ffungene$fungene<-str_split_fixed(ffungene$fungene,';',2)[,1]
  ffungene[2:length(colnames(ffungene))]<-as.data.frame(lapply(ffungene[2:length(colnames(ffungene))],as.numeric))

  raw <-ffungene
  colnames(raw)[1]<-"function"
  processed <- inner_join(raw,KO,by = c("function" = "KO"))

  processed<-subset(processed,select = c(2:ncol(raw),L1,L2,L3,PathwayID))

  processed1<-unite(processed,"com",L1,L2,L3,PathwayID,sep=";")

  x <- processed1%>%
    group_by(com)%>%
    summarise_all(sum)

  x[,c("L1","L2","L3","PathwayID")]<-str_split_fixed(x$com,pattern=";",4)
  x<-x[,-1]

  if (ko1.select=="all") {

    x<-x
    cat("可视化所有层级","\n","\n")
  } else {
    x<-x[which(x$L1==ko1.select),]
    cat("可视化",ko1.select,"层级下功能通路","\n","\n")
  }


  sampleda <- data.frame(sample=colnames(abun),
                         group=substr(colnames(abun),start = 1,stop = 2))%>%column_to_rownames(var = "sample")

  featuretab <- 100000*x[,c(1:(ncol(raw)-1))]
  rownames(featuretab) <- x$PathwayID
  tax <- as.matrix(x[,c(ncol(raw):ncol(x))] %>% tibble::column_to_rownames("PathwayID"))

  library(phyloseq)
  ps <- phyloseq(otu_table(featuretab,taxa_are_rows = T),
                 phyloseq::tax_table(tax),sample_data(sampleda))

  deres <- MicrobiotaProcess::diff_analysis(obj=ps, class="group",
                                            mlfun="lda",
                                            filtermod="pvalue", ##  pvalue
                                            firstcomfun = "kruskal.test",  #  oneway.test
                                            firstalpha=0.05,
                                            strictmod=FALSE,
                                            secondcomfun = "wilcox.test",
                                            subclmin=2,
                                            subclwilc=TRUE,
                                            secondalpha=0.05,
                                            lda=ldascore,
                                            type="others" #非物种
  )

  deres

  p2<-MicrobiotaProcess::ggdiffclade(obj=deres,
                                     layout=layout,#布局类型
                                     alpha=0.5, #树分支的背景透明度
                                     linewd=0.2, #树中连线粗细
                                     skpointsize=0.8, #树骨架中国点的大小
                                     taxlevel=3, #要展示的树分支水平1（L1）及1以下（L2和L3）
                                     cladetext=2,#文本大小
                                     setColors=T,#自定义颜色
                                     removeUnknown=T,#不删分支，但移除分类中有un_的物种注释
                                     reduce=T)+ # 移除分类不明的物种分支，分支和注释一同删去，为T时removeUnknown参数无效
    scale_fill_manual(values=color_group,
                      guide=guide_legend(order=2)) +
    scale_size_continuous(range = c(0.5, 2),
                          guide = guide_legend(keywidth = 0.2,
                                               keyheight = 1,
                                               order = 1))+
    guides(#col = guide_legend(nrow = 8),
      color = guide_legend(keywidth = 0.2, keyheight = 1,order = 2,ncol=1)) +
    theme(panel.background=element_rect(fill=NA),
          legend.position="right",
          plot.margin=margin(0,0,0,0),
          axis.title.x=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = -1),
          axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          legend.spacing.y=unit(0.02, "cm"),
          legend.title=element_text(size=12,face = "bold",family = "serif"),
          legend.text=element_text(size=8,face = "bold",family = "serif"),
          legend.box.spacing=unit(0.02,"cm"))

  ggsave(paste(export_path,"/funpred_cladogram_",ko1.select,".pdf",sep = ""),p2)
  cat("\n","functional_prediction_cladogram.pdf已保存至相对路径","/",export_path,"","下","\n","\n",sep = "")
  return(list(cladogram=p2,pred_res=deres))
}


"tax4Fun"<-function (Tax4FunInput, folderReferenceDatax,
                     fctProfiling = TRUE,
                     refProfile = "UProC",
                     shortReadMode = TRUE,
                     normCopyNo = TRUE)
{
  Tax4FunReferenceData <- microchat::importTax4FunReferenceData(folderReferenceDatax)
  commonOTUs <- intersect(Tax4FunReferenceData$SilvaIDs$V1,
                          row.names(Tax4FunInput$otuTable))
  indexInput <- match(commonOTUs, row.names(Tax4FunInput$otuTable))
  indexSILVAToKEGG <- match(commonOTUs, Tax4FunReferenceData$SilvaIDs$V1)
  subsetOTUTables <- as.matrix(Tax4FunInput$otuTable[indexInput,
  ])
  subsetSILVAToKEGG <- Tax4FunReferenceData$SilvaToKEGGMappingMat[indexSILVAToKEGG,
  ]
  subsetSILVAToKEGG <- as.data.frame(as.matrix(subsetSILVAToKEGG))
  if (fctProfiling) {
    FctCat <- Tax4FunReferenceData$KEGGKOInformation
    if (refProfile == "UProC") {
      if (shortReadMode) {
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCShort
      } else {
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCLong
      }
    } else if (refProfile == "PAUDA") {
      if (shortReadMode) {
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchPAUDAShort
      } else {
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchPAUDALong
      }
    } else {
      print("Invalid functional profiling method. Using default functional profiling method.")
      if (shortReadMode) {
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCShort
      } else {
        RefProfile <- Tax4FunReferenceData$FctAbundancesKEGGBacArchUProCLong
      }
    }
  } else {
    RefProfile <- Tax4FunReferenceData$PathwayAbundancesKEGGBacArch
    FctCat <- Tax4FunReferenceData$PathwayInformationKEGGBacArch
  }
  Tax4FunProfile <- matrix(data = 0, nrow = ncol(subsetOTUTables),
                           ncol = nrow(RefProfile))
  for (i in 1:ncol(subsetOTUTables)) {
    KEGGTaxProfile <- t(subsetSILVAToKEGG) %*% subsetOTUTables[,i]
    if (normCopyNo) {
      NormKEGGTaxProfile <- KEGGTaxProfile/Tax4FunReferenceData$KEGGBacArchCopyNumbers
    } else {
      NormKEGGTaxProfile <- as.data.frame(KEGGTaxProfile)
    }
    NormKEGGTaxProfile <- NormKEGGTaxProfile/sum(NormKEGGTaxProfile)
    Tax4FunProfile[i, ] <- as.matrix(RefProfile) %*% NormKEGGTaxProfile$V1
  }
  if (fctProfiling) {
    FctCat <- FctCat[which(!colSums(Tax4FunProfile) == 0),
                     c(1, 2)]
    Tax4FunProfile <- Tax4FunProfile[, which(!colSums(Tax4FunProfile) ==
                                               0)]
    Tax4FunProfile <- as.matrix(Tax4FunProfile)
    if (ncol(Tax4FunProfile) > 1) {
      colnames(Tax4FunProfile) <- paste(FctCat$V2, FctCat$V1,
                                        sep = "; ")
      rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
    } else {
      Tax4FunProfile <- as.matrix(t(Tax4FunProfile))
      colnames(Tax4FunProfile) <- paste(FctCat$V2, FctCat$V1,
                                        sep = "; ")
      rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
    }
  } else {
    FctCat <- FctCat[which(!colSums(Tax4FunProfile) == 0), c(1, 2)]
    Tax4FunProfile <- Tax4FunProfile[, which(!colSums(Tax4FunProfile) == 0)]
    Tax4FunProfile <- as.matrix(Tax4FunProfile)
    if (ncol(Tax4FunProfile) > 1) {
      colnames(Tax4FunProfile) <- paste(FctCat$V3, FctCat$V4,
                                        sep = "; ")
      rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
    } else {
      Tax4FunProfile <- as.matrix(t(Tax4FunProfile))
      colnames(Tax4FunProfile) <- paste(FctCat$V3, FctCat$V4,
                                        sep = "; ")
      rownames(Tax4FunProfile) <- Tax4FunInput$sampleNames
    }
  }
  FTU <- 1 - colSums(subsetOTUTables)/colSums(Tax4FunInput$otuTable)
  names(FTU) <- Tax4FunInput$sampleNames
  Tax4FunProfile <- list(Tax4FunProfile = Tax4FunProfile,
                         FTU = FTU, fctProfiling = fctProfiling, refProfile = refProfile,
                         shortReadMode = shortReadMode)
  class(Tax4FunProfile) <- "microchat"
  return(Tax4FunProfile)
}


"importTax4FunReferenceData"<-function (folder)
{
  if (substr(folder, nchar(folder), nchar(folder)) == "/") {
    pathReferenceData <- folder
  }
  else {
    pathReferenceData <- paste(folder, "/", sep = "")
  }
  referenceData <- list()
  tmpReferenceData <- readRDS(paste(pathReferenceData, "KEGGBacArchTaxInformationMoPPro.RData",
                                    sep = ""))
  referenceData$KEGGBacArchTaxInformationMoPPro <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "PathwayAbundancesKEGGBacArchMoPPro.RData",
                                    sep = ""))
  referenceData$PathwayAbundancesKEGGBacArchMoPPro <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "PathwayInformationKEGGBacArchMoPPro.RData",
                                    sep = ""))
  referenceData$PathwayInformationKEGGBacArchMoPPro <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "PathwayAbundancesKEGGBacArch.RData",
                                    sep = ""))
  referenceData$PathwayAbundancesKEGGBacArch <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "PathwayInformationKEGGBacArch.RData",
                                    sep = ""))
  referenceData$PathwayInformationKEGGBacArch <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "KEGGKOInformation.RData",
                                    sep = ""))
  referenceData$KEGGKOInformation <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "KEGGBacArchTaxInformation.RData",
                                    sep = ""))
  referenceData$KEGGBacArchTaxInformation <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "KEGGBacArchCopyNumbers.RData",
                                    sep = ""))
  referenceData$KEGGBacArchCopyNumbers <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "FctAbundancesKEGGBacArchPAUDALong.RData",
                                    sep = ""))
  referenceData$FctAbundancesKEGGBacArchPAUDALong <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "FctAbundancesKEGGBacArchPAUDAShort.RData",
                                    sep = ""))
  referenceData$FctAbundancesKEGGBacArchPAUDAShort <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "FctAbundancesKEGGBacArchUProCLong.RData",
                                    sep = ""))
  referenceData$FctAbundancesKEGGBacArchUProCLong <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "FctAbundancesKEGGBacArchUProCShort.RData",
                                    sep = ""))
  referenceData$FctAbundancesKEGGBacArchUProCShort <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "SilvaToKEGGMappingMat.RData",
                                    sep = ""))
  referenceData$SilvaToKEGGMappingMat <- tmpReferenceData
  tmpReferenceData <- readRDS(paste(pathReferenceData, "SilvaIDs.RData",
                                    sep = ""))
  referenceData$SilvaIDs <- tmpReferenceData
  return(referenceData)
}


"calcTaxFun"<-function(submchat,
                       comparison_sel="ct-xl",
                       silva_path='/Volumes/BOOTCAMP/SILVA123',
                       export_path="tax4fun") {
  export_path<-paste(export_path,"/microbial functional prediction analysis/data/",comparison_sel,sep = "")
  dir.create(export_path, recursive = TRUE)

  data(kegg_anno)
  kegg_anno=kegg_anno

  qiimeotu <- readtaxdata(submchat,comparison = comparison_sel)
  QIIMESingleData<-tax4funpretreat(qiimeotu)
  tax4fp<-tax4funpred(QIIMESingleData=QIIMESingleData,
                      kegg_anno=kegg_anno,
                      folderReferenceData = silva_path,
                      export_path=export_path)
  group<-qiimeotu$sampleda
  fungene<-tax4fp$tax4fun_gene%>%rownames_to_column(var = "fungene")
  funpred<-tax4fp$tax4fun_pathway

  return(list(funpred=funpred,
              fungene=fungene,
              qiimeotu=qiimeotu,
              submchat=submchat,
              comparison_sel=comparison_sel))
}

"calcTaxFunDiff.single" <- function(ko1.select="Metabolism", ###  ko1.select="all",则使用所有ko1
                                    microchatTaxFun,
                                    padj=FALSE,
                                    min.fun.abun = 0.3,   ###平均丰度大于该值才考虑组图
                                    filter=FALSE,          ###是否过滤掉差异不显著的
                                    export_path) {

  fp.dat<-microchatTaxFun
  abun<-fp.dat$submchat$otu_table
  tt<-str_split_fixed(fp.dat$comparison_sel,'-',2)%>%as.character()
  message("Control selected is ",tt[1])
  message("Treatment selected is ",tt[2])
  export_path<-paste(export_path,"/microbial functional prediction analysis/Pathway/",tt[1],"-",tt[2],sep = "")
  dir.create(export_path, recursive = TRUE)

  funpred <- fp.dat$funpred
  qiimeotu <- fp.dat$qiimeotu
  group<-qiimeotu$sampleda

  a<-subset(funpred,select=c(2:(nrow(group)+1),KO1,KO2,KO3))
  datako2<-subset(funpred,select=c(2:(nrow(group)+1),KO2))
  datako3<-subset(funpred,select=c(2:(nrow(group)+1),KO3))
  datako4<-subset(funpred,select=c(2:(nrow(group)+1),KO1,KO3))
  datako4<-unite(datako4,"KO3",KO1,KO3,sep = ";")

  datako4 %>%
    group_by(KO3) %>%
    summarise_all(sum) -> data
  data[,c("KO1","KO3")]<-str_split_fixed(data$KO3,pattern=";",2)
  data<-data[which(data$KO1==ko1.select),]

  data<-subset(data,select = c(-KO1))
  data<-column_to_rownames(data,var = "KO3")

  data <- data*100
  data <- data %>% filter(apply(data,1,mean) > min.fun.abun)
  data <- t(data)
  data1 <- data.frame(data,group$group)
  colnames(data1) <- c(colnames(data),"Group")
  data1$Group <- as.factor(data1$Group)



  diff <- data1 %>%
    select_if(is.numeric) %>%
    map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')

  if (padj) diff$p.value <- p.adjust(diff$p.value,"fdr")
  if (filter) {
    diff1 <- diff %>% filter(p.value < 0.05)
    if (nrow(diff1) <= 0) {
      diff<-diff
    } else {
      diff<-diff1
    }
  }

  ## left barplot layout
  abun.bar <- data1[,c(diff$var,"Group")] %>%
    gather(variable,value,-Group) %>%
    group_by(variable,Group) %>%
    summarise(Mean = mean(value))

  kk1<-subset(funpred,select=c(KO1,KO2,KO3))
  kk1$variable<-kk1$KO3
  KKK<-merge(abun.bar,kk1,by="variable")
  abun.bar<-KKK

  ## right point layout
  diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
  diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                              levels(data1$Group)[2]))
  diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]

  abun.bar<-abun.bar[order(abun.bar$KO1,abun.bar$KO2),]
  #abun.bar<-abun.bar[order(abun.bar$KO2),]
  rownames(abun.bar)<-1:nrow(abun.bar)
  abun.bar$variable<-factor(abun.bar$variable,levels=rev(unique(abun.bar$variable)))

  diff.mean$KO3<-diff.mean$var
  ajj<-merge(diff.mean,kk1,by="KO3", all = TRUE,sort = FALSE)
  df<-filter(ajj,ajj$var!="NA")
  abun.bar1<-abun.bar

  abun.bar1$KO4<-abun.bar1$KO2
  abun.bar1$KO4[setdiff(1:nrow(abun.bar1),match(unique(abun.bar1$KO4),abun.bar1$KO4))]<-""
  abun.bar1$KO4[match(unique(abun.bar1$KO4),abun.bar1$KO4)]<-abun.bar1$KO4[match(unique(abun.bar1$KO4),abun.bar1$KO4)]


  ## right text layout
  abun.bar33<-subset(abun.bar1,select=c(KO2,KO3))
  abun.bar33<-abun.bar33[duplicated(abun.bar33$KO3),]
  rownames(abun.bar33)<-1:nrow(abun.bar33)

  diff.mean1<-diff.mean
  diff.mean1[,"KO2"]<-abun.bar33$KO2[match(diff.mean$KO3,abun.bar33$KO3)]
  diff.mean1<-diff.mean1[order(diff.mean1$KO2,diff.mean1$KO3),]

  diff.mean1$KO2 <- factor(diff.mean1$KO2,levels = rev(unique(diff.mean1$KO2)))
  diff.mean1$p.value <- signif(diff.mean1$p.value,3)
  diff.mean1$p.value <- as.numeric(diff.mean1$p.value)
  diff.mean1$var<-factor(diff.mean1$var,levels =rev(diff.mean1$var))

  diff.mean1[,"sig"]<-ifelse(diff.mean1$p.value<0.001,"***",
                             ifelse(diff.mean1$p.value<0.01,"**",
                                    ifelse(diff.mean1$p.value<0.05,"*","ns") ))

  cat("Totally seeking ",nrow(funpred)," (KO3) pathways from ",length(unique(funpred$KO2))," (KO2) pathways.\n",sep = "")
  cat("Totally seeking ",nrow(diff.mean1)," (KO3) significant pathways, accounted for ",paste(sprintf("%0.2f",round(100*nrow(diff.mean1)/nrow(funpred))),"%",sep = ""),".\n",sep = "")
  cat("Totally seeking ",length(unique(diff.mean1$KO2))," (KO2) significant pathways, accounted for ",paste(sprintf("%0.2f",round(100*length(unique(diff.mean1$KO2))/length(unique(funpred$KO2)))),"%",sep = ""),".\n",sep = "")

  file2=paste(export_path, "/tax4fun_FunPathway",ko1.select,"_",paste(tt[1],tt[2],sep = "-"),"_mean.txt",sep = "" )
  write.table(abun.bar,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Differential analysis results based on welch's t test using tax4fun prediction (",paste(tt[1],tt[2],sep = "-"),") has been exported to ","/",export_path,"",sep = "","\n")

  file2=paste(export_path, "/tax4fun_FunPathway",ko1.select,"_",paste(tt[1],tt[2],sep = "-"),"_onlydiffdata.txt",sep = "" )
  write.table(diff.mean1,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Differential analysis results based on welch's t test using tax4fun prediction (",paste(tt[1],tt[2],sep = "-"),") has been exported to ","/",export_path,"",sep = "","\n")

  return(list(abun.bar=abun.bar,diff.mean=diff.mean,df=df,
              comparison=fp.dat$comparison_sel,otu_table=abun,
              abun.bar1=abun.bar1,diff.mean1=diff.mean1,sampledata=microchatTaxFun$qiimeotu$sampleda))
}

"calcTaxFunDiff.all" <- function(microchatTaxFun,
                                 padj=FALSE,
                                 min.fun.abun,
                                 filter=FALSE,
                                 export_path=export_path) {
  ddxdata<-lapply(unique(microchatTaxFun$funpred$KO1)[c(6,4:1)],
                  FUN = calcTaxFunDiff,
                  microchatTaxFun=microchatTaxFun,
                  padj=padj,
                  min.fun.abun = min.fun.abun,
                  filter=filter,
                  export_path=export_path)

  abun.barx<-data.frame()
  diff.meanx<-data.frame()
  dfx<-data.frame()
  abun.bar1x<-data.frame()
  diff.mean1x<-data.frame()
  for (i in 1:length(ddxdata)) {
    dat.sel<-ddxdata[[i]]

    abun.bar<-dat.sel$abun.bar
    diff.mean<-dat.sel$diff.mean
    df<-dat.sel$df
    abun.bar1<-dat.sel$abun.bar1
    diff.mean1<-dat.sel$diff.mean1
    {
      abun.barx<-rbind(abun.barx,abun.bar)
      diff.meanx<-rbind(diff.meanx,diff.mean)
      dfx<-rbind(dfx,df)
      abun.bar1x<-rbind(abun.bar1x,abun.bar1)
      diff.mean1x<-rbind(diff.mean1x,diff.mean1)
    }
  }

  microchatTaxFun.all<-list(abun.bar=abun.barx,
                            diff.mean=diff.meanx,
                            df=dfx,
                            abun.bar1=abun.bar1x,
                            diff.mean1=diff.mean1x,
                            otu_table=ddxdata[[1]]$otu_table,
                            sampledata=ddxdata[[1]]$sampledata,
                            comparison=ddxdata[[1]]$comparison)
  return(microchatTaxFun.all)
}

"calcTaxFunDiff"<-function(ko1.select="Metabolism", ###  ko1.select="all",则使用所有ko1
                           microchatTaxFun,
                           padj=FALSE,
                           min.fun.abun = 0.3,   ###平均丰度大于该值才考虑组图
                           filter=FALSE,          ###是否过滤掉差异不显著的
                           export_path) {
  if (ko1.select!="all") {
    microtfdiff<-calcTaxFunDiff.single(
      ko1.select=ko1.select,
      microchatTaxFun=microchatTaxFun,
      padj=padj,
      min.fun.abun = min.fun.abun,
      filter=filter,
      export_path=export_path)
    cat("\n","Export KO1 pathway: ",ko1.select)
  } else {
    microtfdiff<-calcTaxFunDiff.all(microchatTaxFun=microchatTaxFun,
                                    padj=padj,
                                    min.fun.abun=min.fun.abun,
                                    filter=filter,
                                    export_path=export_path)
    cat("\n","Export all KO1 pathway.")
  }

  return(microtfdiff)
}

"gradp" <- function(color.sel,abun.bar,p1,abun.bar1,ko2.bar.alpha) {
  colorsbar<-color.sel
  abun.barzz<-abun.bar
  ko2.selected<-c()
  for (i in unique(abun.barzz$KO1)) {
    ko.selected<-abun.barzz[which(abun.barzz$KO1==i),]$KO2%>%unique()%>%rev()
    {
      ko2.selected<-c(ko2.selected,ko.selected)
    }
  }
  suppressMessages(ydensx <- cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
                     ggnewscale::new_scale_fill()+
                     ggplot2::scale_x_discrete()+
                     ggplot2::coord_flip()+
                     ggplot2::geom_col(data=abun.bar1, ggplot2::aes(x=variable,
                                                                    y=1,
                                                                    fill=KO2),
                                       color=NA,
                                       alpha=ko2.bar.alpha, width = 1,
                                       size=.2) +
                     ggplot2::scale_fill_manual(values = colorsbar)+
                     ggplot2::theme(panel.spacing = ggplot2::unit(0,"lines"))+
                     ggplot2::labs(title="KO2   KO3"))
  return(ydensx)
}




"plotDiffPathwayExtendbar" <- function(microchatTaxFunDiff,
                                       xlabname=xlabname,
                                       ko2.bar.gradient=FALSE,
                                       ko2.bar.gradient.rev=FALSE,
                                       gradient.num=20,
                                       sig.text.size=3,
                                       sig.text.color="black",
                                       ko2.text.size=2,
                                       ko2.text.color="black",
                                       add.bar.text.size=3,
                                       add.bar.text.color="white",
                                       pdf.width=NULL,
                                       layout.rt=c(1,4,4,4,1.5),
                                       ko1.color=colorCustom(10,pal = "set1"),
                                       add.bar.color=colorCustom(50,pal = "gygn"),
                                       compare.color=colorCustom(4,pal = "gygn"),
                                       ko1.bar.alpha=1,
                                       ko2.bar.alpha=1,
                                       ko1.text.size=4,
                                       ko1.text.color="black",
                                       siglabel=c("sig","label","sig+label"),
                                       export_path=export_path) {

  siglabel<-match.arg(siglabel)
  abun.bar<-abun.bar<-microchatTaxFunDiff$abun.bar
  diff.mean<-microchatTaxFunDiff$diff.mean
  df<-microchatTaxFunDiff$df
  abun.bar1<-microchatTaxFunDiff$abun.bar1
  diff.mean1<-microchatTaxFunDiff$diff.mean1
  sampledata<-microchatTaxFunDiff$sampledata
  comparison=microchatTaxFunDiff$comparison
  abun<-microchatTaxFunDiff$otu_table


  export_path<-paste(export_path,"/microbial functional prediction analysis",sep = "")
  dir.create(export_path, recursive = TRUE)

  tt<-str_split_fixed(comparison,'-',2)%>%as.character()
  color.all<-compare.color
  allg<-unique(substr(colnames(abun),start = 1,stop = 2))
  if (length(color.all)<length(allg)) stop("Please provide more group colors, which is ",length(allg),".")
  color.use<-color.all[match(tt,allg)]

  gp.use<-xlabname[match(tt,allg)]

  cbbPalette <- color.use
  names(cbbPalette)<-unique(sampledata$group)

  p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
    scale_x_discrete(limits = levels(abun.bar$variable),label=NULL) +
    scale_fill_manual(values=cbbPalette,
                      labels=gp.use,
                      limits = (unique(abun.bar$Group)))+
    coord_flip() +
    xlab("") +
    ylab("Mean proportion (%)") +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length.y = unit(0,"lines"),
          axis.ticks.y  = element_line(size=0),
          axis.line  = element_line(colour = "black"),
          axis.title.y=element_text(colour='black', size=8,face = "bold"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          legend.title=element_blank(),
          legend.background = element_blank(),
          legend.text=element_text(size=8,face = "bold",colour = "black",
                                   family = "serif",
                                   margin = margin(r = 20)),
          legend.position = c(-0.5,-0.025),
          legend.direction = "horizontal",
          legend.key = element_rect(colour = "white"),
          legend.key.width = unit(0.8,"cm"),
          legend.key.height = unit(0.5,"cm"))

  for (i in 1:(nrow(diff.mean) - 1)) {
    p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }


  p1 <- p1 +
    geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black")

  if (length(add.bar.color)<length(unique(abun.bar$KO2))) stop("Please provide more colors for coloring KEGG 2 levels, which is ",length(unique(abun.bar$KO2)),".")
  colorsbar<-add.bar.color[1:length(unique(abun.bar$KO2))]

  abun.barzz<-abun.bar
  ko2.selected<-c()
  for (i in unique(abun.barzz$KO1)) {
    ko.selected<-abun.barzz[which(abun.barzz$KO1==i),]$KO2%>%unique()%>%rev()
    {
      ko2.selected<-c(ko2.selected,ko.selected)
    }
  }
  names(colorsbar)<-ko2.selected%>%rev()
  # need to set `coord_flip = TRUE` if you plan to use `coord_flip()`
  if (!ko2.bar.gradient) {
    ydens <- cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
      ggnewscale::new_scale_fill()+
      scale_x_discrete()+
      coord_flip()+
      geom_col(data=abun.bar1, aes(x=variable,
                                   y=1,
                                   fill=KO2),
               color=NA,
               alpha=ko2.bar.alpha, width = 1,
               size=.2) +
      geom_text(data=abun.bar1, aes(x=variable,
                                    y=2,
                                    label=variable),size=add.bar.text.size,
                color=add.bar.text.color,
                family="serif",check_overlap = FALSE,
                hjust = 1)+
      geom_text(data=abun.bar1, aes(x=variable,
                                    y=0,
                                    label=KO4),size=ko2.text.size,
                color=ko2.text.color,
                family="serif",check_overlap = FALSE,
                hjust = 0)+
      scale_fill_manual(values = colorsbar)+
      labs(title="KO2   KO3")
  } else {
    colorsbar<-add.bar.color[1:length(unique(abun.bar$KO2))]

    colorss<-list()
    for (j in 1:length(colorsbar)) {
         if (!ko2.bar.gradient.rev) {
           colorss[[j]] <- grDevices::colorRampPalette(c("white",
                                                         colorsbar[j]))(gradient.num)
         } else {
           colorss[[j]] <- grDevices::colorRampPalette(c(colorsbar[j],
                                                         "white"))(gradient.num)
      }

    }

    color.sel<-list()
    for (j in 1:length(colorss[[1]])) {
      c.sel<-c()
      for (i in 1:length(colorss)) {
        c.sel<-c(c.sel,colorss[[i]][j])
        color.sel[[j]]<-c.sel
      }
    }

    nworker<-parallel::detectCores()-2
    if (nworker<=0) nworker<-parallel::detectCores()
    if (nworker!=1) {
      c1 <- try(parallel::makeCluster(nworker, type = "PSOCK"))
      if (class(c1)[1] == "try-error") {
        c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))
      }
      if (class(c1)[1] == "try-error") {
        c1 <- try(parallel::makeCluster(nworker, setup_strategy = "sequential"))
      }
      message("\n","Now generating gradient color by multi threaded computing. Begin at ",
              date(), ". Please wait...")
      clusterExport(c1, "%>%")
      ydensx = parallel::parLapply(c1, color.sel, gradp,
                                   abun.bar,p1,abun.bar1,ko2.bar.alpha)
      parallel::stopCluster(c1)
    }

    ydensx<-patchwork::wrap_plots(ydensx,nrow = 1,byrow = TRUE)

    ydens2 <- cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
      ggnewscale::new_scale_fill()+
      scale_x_discrete()+
      coord_flip()+
      geom_col(data=abun.bar1, aes(x=variable,
                                   y=1), alpha=ko2.bar.alpha, width = 1,
               fill=NA,
               size=.2) +
      geom_text(data=abun.bar1, aes(x=variable,
                                    y=0,
                                    label=KO4),size=ko2.text.size,
                color=ko2.text.color,inherit.aes = TRUE,
                family="serif",check_overlap = FALSE,
                hjust = 0)+
      labs(title="KO2   KO3")

    ydens3 <- cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
      ggnewscale::new_scale_fill()+
      scale_x_discrete()+
      coord_flip()+
      geom_col(data=abun.bar1, aes(x=variable,
                                   y=1),
               fill=NA,
               alpha=ko2.bar.alpha, width = 1,
               size=.2) +
      geom_text(data=abun.bar1, aes(x=variable,
                                    y=2,
                                    label=variable),
                size=add.bar.text.size,
                color=add.bar.text.color,inherit.aes = TRUE,
                family="serif",check_overlap = FALSE,
                hjust = 1)+
      labs(title="KO2   KO3")

    ydensxZ<-ydensx+annotation_custom(grob=ggplotGrob(ydens3),
                                      ymin = -(gradient.num*2)+2, ymax = 2,
                                      xmin=0.4,
                                      xmax=0.6+nrow(abun.bar)/2)

    ydens<-ydensxZ+annotation_custom(grob=ggplotGrob(ydens2),
                                     ymin = -(gradient.num*2)+2, ymax = 2,
                                     xmin=0.4,
                                     xmax=0.6+nrow(abun.bar)/2)

  }

  p2 <- ggplot(diff.mean1,aes(x=var,y=estimate,fill = Group)) +
    scale_x_discrete(limits = levels(abun.bar$variable)) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +

    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals")


  for (i in 1:(nrow(diff.mean1) - 1))
    p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

  p2 <- p2 +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),show.legend = FALSE,
                  position = position_dodge(0.8), width = 0.5, size = 0.5) +
    geom_point(shape = 21,size = 3) +
    scale_fill_manual(values=cbbPalette) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')



  if (siglabel=="label"){
    diff.mean1$p.value<-diff.mean1$p.value%>%as.numeric()
    diff.mean1$p.value<-sprintf("%0.3f",round(diff.mean1$p.value,3))
    p3 <- ggplot(diff.mean1,aes(var,estimate,fill = Group)) +
      scale_x_discrete(limits = levels(abun.bar$variable)) +
      geom_text(aes(y = 0,x = var),label = diff.mean1$p.value,family="serif",
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = sig.text.size) +
      geom_text(aes(x = nrow(diff.mean1)/2 +0.5,y = 0.85),label = "P-value (corrected)",
                family="serif",srt = 90,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_text(family = "serif"),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }

  if (siglabel=="sig"){
    diff.mean1$p.value<-diff.mean1$p.value%>%as.numeric()
    diff.mean1$p.value<-sprintf("%0.3f",round(diff.mean1$p.value,3))

    p3 <- ggplot(diff.mean1,aes(var,estimate,fill = Group)) +
      scale_x_discrete(limits = levels(abun.bar$variable)) +
      geom_text(aes(y = 0,x = var),label = diff.mean1$sig,family="serif",
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = sig.text.size) +
      geom_text(aes(x = nrow(diff.mean1)/2 +0.5,y = 0.85),label = "P-value (corrected)",
                family="serif",srt = 90,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_text(family = "serif"),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }

  if (siglabel=="sig+label"){
    diff.mean1$p.value<-diff.mean1$p.value%>%as.numeric()
    diff.mean1$p.value<-sprintf("%0.3f",round(diff.mean1$p.value,3))

    p3 <- ggplot(diff.mean1,aes(var,estimate,fill = Group)) +
      scale_x_discrete(limits = levels(abun.bar$variable)) +
      geom_text(aes(y = 0,x = var),label = paste(diff.mean1$p.value,diff.mean1$sig,sep = ""),
                family="serif",color=sig.text.color,show.legend = FALSE,
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = sig.text.size) +
      geom_text(aes(x = nrow(diff.mean1)/2 +0.5,y = 0.85),
                label = "P-value (corrected)",show.legend = FALSE,
                family="serif",srt = 270,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_text(family = "serif"),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }

  ## merge plots
  library(patchwork)
  ko1.color<-ko1.color[1:length(unique(abun.bar1$KO1))]
  names(ko1.color)<-unique(abun.bar1$KO1)%>%rev()

  ko1.text<-table(abun.bar1$KO1)
  ko1.text<-ko1.text[match(names(ko1.color),names(ko1.text))]
  ko1.text<-rev(ko1.text)
  ko1.text<-ko1.text%>%data.frame()
  if (nrow(ko1.text)==1) {
    ko1.text$Freq<-ko1.text$.
    ko1.text$freq<-ko1.text$Freq

    ko1.text$half<-ko1.text$Freq/2
    ko1.text$ycoord<-ko1.text$freq-ko1.text$half

    abun.bar1$KO5<-""
    ypos<-ko1.text$ycoord
    abun.bar1$KO5[ypos]<-rownames(ko1.text)%>%as.character()
  } else {
    ko1.text$freq<-ko1.text$Freq%>%cumsum()
    ko1.text$half<-ko1.text$Freq/2
    ko1.text$ycoord<-ko1.text$freq-ko1.text$half

    abun.bar1$KO5<-""
    for (i in 1:nrow(ko1.text)) {
      ypos<-ko1.text$ycoord[i]
      abun.bar1$KO5[ypos]<-ko1.text$Var1[i]%>%as.character()
    }
  }

  #str_split_fixed("Human Diseases"," ",n=2)
  abun.bar1$KO5[which(abun.bar1$KO5=="Human Diseases")]<-"HD"
  abun.bar1$KO5[which(abun.bar1$KO5=="Genetic Information Processing")]<-"GIP"
  abun.bar1$KO5[which(abun.bar1$KO5=="Environmental Information Processing")]<-"EIP"
  abun.bar1$KO5[which(abun.bar1$KO5=="Cellular Processes")]<-"CellPro"

  p4 <- cowplot::axis_canvas(p1, axis = "y", coord_flip = TRUE) +
    ggnewscale::new_scale_fill()+
    scale_x_discrete()+
    coord_flip()+
    geom_col(data=abun.bar1, aes(x=variable,
                                 y=1,
                                 fill=KO1),
             alpha=ko1.bar.alpha,
             width = 1,show.legend = FALSE,
             size=.2) +
    scale_fill_manual(values = ko1.color)


  whatsn<-abun.bar1$KO5[which(abun.bar1$KO5!="")]
  whats<-(nrow(abun.bar1)/2)-(which(abun.bar1$KO5!="")/2)
  for (i in 1:length(whats)) {
    p4<-p4+annotate(geom = "text",
                    family="serif",
                    srt = 90,
                    size=ko1.text.size,
                    color=ko1.text.color,
                    x = whats[i]+0.5,
                    y = 1,
                    label = whatsn[i])
  }

  #p4+ydens+plot_layout(widths = c(0.5,0.5),byrow=TRUE)
  if (is.null(pdf.width)) widthx=15 else widthx= pdf.width
  heighty=length(unique(abun.bar1$KO3))/4
  if (heighty<(widthx/2)) heighty<-widthx/2 else  heighty<-heighty

  p<-p4+ydens+p1 + p2 + p3 + plot_layout(widths = c(layout.rt[1],
                                                    layout.rt[2],
                                                    layout.rt[3],
                                                    layout.rt[4],
                                                    layout.rt[5]),byrow=TRUE)

  message("width used for ggsave is ",widthx,", as well as ",heighty," of height")


  ggsave(paste(export_path,"/(",comparison,") Extended Errorbar plot (CMYK).pdf",sep = ""),
         colormodel="cmyk",
         width = widthx,height = heighty,p)
  ggsave(paste(export_path,"/(",comparison,") Extended Errorbar plot (RGB).pdf",sep = ""),
         colormodel="srgb",
         width = widthx,height = heighty,p)
  ggsave(paste(export_path,"/(",comparison,") Extended Errorbar plot.tiff",sep = ""),
         bg="white",width = widthx,height = heighty,p)
  cat("\n已输出",
      paste(comparison,"基于","welch's t test","的Extended error bar plot",sep = ""))
  print(p)
  return(p)
}

"plotDiffGeneExtendbar" <- function(microchatTaxFunGeneDiff,
                                    xlabname=xlabname,
                                    sig.text.size=3,
                                    sig.text.color="black",
                                    pdf.width=NULL,
                                    layout.rt=c(4,4,4),
                                    add.bar.color=colorCustom(50,pal = "gygn"),
                                    compare.color=colorCustom(4,pal = "gygn"),
                                    export_path=export_path) {

  abun.bar<-abun.bar<-microchatTaxFunGeneDiff$abun.bar
  diff.mean<-microchatTaxFunGeneDiff$diff.mean
  diff.mean1<-microchatTaxFunGeneDiff$diff.mean1
  sampledata<-microchatTaxFunGeneDiff$sampledata
  comparison=microchatTaxFunGeneDiff$comparison
  abun<-microchatTaxFunGeneDiff$otu_table

  export_path<-paste(export_path,"/microbial functional prediction analysis/Functional gene/",sep = "")
  dir.create(export_path, recursive = TRUE)

  tt<-str_split_fixed(comparison,'-',2)%>%as.character()
  color.all<-compare.color
  allg<-unique(substr(colnames(abun),start = 1,stop = 2))
  if (length(color.all)<length(allg)) stop("Please provide more group colors, which is ",length(allg),".")
  color.use<-color.all[match(tt,allg)]

  gp.use<-xlabname[match(tt,allg)]

  cbbPalette <- color.use
  names(cbbPalette)<-unique(sampledata$group)
  abun.bar$variable
  abun.bar[,c("var","name")]<-str_split_fixed(abun.bar$variable,";",2)
  abun.bar$var<-factor(abun.bar$var,levels = unique(abun.bar$var))
  p1 <- ggplot(abun.bar,aes(var,Mean,fill = Group)) +
    scale_x_discrete(limits = levels(abun.bar$var)) +
    scale_fill_manual(values=cbbPalette,
                      labels=gp.use,
                      limits = rev(unique(abun.bar$Group)))+
    coord_flip() +
    xlab("") +
    ylab("Mean proportion (%)") +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length.y = unit(0,"lines"),
          axis.ticks.y  = element_line(size=0),
          axis.line  = element_line(colour = "black"),
          axis.title.y=element_text(colour='black', size=8,face = "bold"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          legend.title=element_blank(),
          legend.background = element_blank(),
          legend.text=element_text(size=8,face = "bold",colour = "black",
                                   family = "serif",
                                   margin = margin(r = 20)),
          legend.position = c(0.2,1.01),
          legend.direction = "horizontal",
          legend.key = element_rect(colour = "white"),
          legend.key.width = unit(0.8,"cm"),
          legend.key.height = unit(0.5,"cm"))

  for (i in 1:(nrow(diff.mean) - 1)) {
    p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }


  p1 <- p1 +
    geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black")

  if (length(add.bar.color)<length(unique(abun.bar$KO2))) stop("Please provide more colors for coloring KEGG 2 levels, which is ",length(unique(abun.bar$KO2)),".")
  colorsbar<-add.bar.color[1:length(unique(abun.bar$KO2))]



  p2 <- ggplot(diff.mean1,aes(x=var,y=estimate,fill = Group)) +
    scale_x_discrete(limits = levels(diff.mean1$var)) +
    coord_flip() +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +

    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals")


  for (i in 1:(nrow(diff.mean1) - 1))
    p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

  p2 <- p2 +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),show.legend = FALSE,
                  position = position_dodge(0.8), width = 0.5, size = 0.5) +
    geom_point(shape = 21,size = 3) +
    scale_fill_manual(values=cbbPalette) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')


  diff.mean1$p.value<-diff.mean1$p.value%>%as.numeric()
  diff.mean1$p.value<-sprintf("%0.3f",round(diff.mean1$p.value,3))

  p3 <- ggplot(diff.mean1,aes(var,estimate,fill = Group)) +
    geom_text(aes(y = 0,x = var),label = paste(diff.mean1$p.value,diff.mean1$sig,sep = ""),
              family="serif",color=sig.text.color,show.legend = FALSE,
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = sig.text.size) +
    geom_text(aes(x = nrow(diff.mean1)/2 +0.5,y = 0.85),
              label = "P-value (corrected)",show.legend = FALSE,
              family="serif",srt = 270,fontface = "bold",size = 5) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          text = element_text(family = "serif"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())

  ## merge plots
  library(patchwork)

  if (is.null(pdf.width)) widthx=15 else widthx= pdf.width
  heighty=length(unique(abun.bar$fungene))/4
  if (heighty<(widthx/2)) heighty<-widthx/2 else  heighty<-heighty

  p<-p1 + p2 + p3 + plot_layout(widths = c(
    layout.rt[1],
    layout.rt[2],
    layout.rt[3]),byrow=TRUE)

  message("width used for ggsave is ",widthx,", as well as ",heighty," of height")

  ggsave(paste(export_path,"/(",comparison,") Fungene Extended Errorbar plot (CMYK).pdf",sep = ""),
         colormodel="cmyk",
         width = widthx,height = heighty,p)
  ggsave(paste(export_path,"/(",comparison,") Fungene Extended Errorbar plot (RGB).pdf",sep = ""),
         colormodel="srgb",
         width = widthx,height = heighty,p)
  ggsave(paste(export_path,"/(",comparison,") Fungene Extended Errorbar plot.tiff",sep = ""),
         width = widthx,height = heighty,p)

  cat("\n已输出",
      paste(comparison,"基于","welch's t test","的Extended error bar plot",sep = ""))
  print(p)
}

"calcTaxFunGeneDiff" <- function(microchatTaxFun,
                                 padj=FALSE,
                                 min.fun.abun = 0.1,   ###平均丰度大于该值才考虑组图
                                 filter=FALSE,          ###是否过滤掉差异不显著的
                                 export_path) {

  fp.dat<-microchatTaxFun
  abun<-fp.dat$submchat$otu_table
  tt<-str_split_fixed(fp.dat$comparison_sel,'-',2)%>%as.character()
  message("Control selected is ",tt[1])
  message("Treatment selected is ",tt[2])
  export_path<-paste(export_path,"/microbial functional prediction analysis/Functional gene/",tt[1],"-",tt[2],sep = "")
  dir.create(export_path, recursive = TRUE)
  qiimeotu<-microchatTaxFun$qiimeotu
  group<-qiimeotu$sampleda
  fungene<-fp.dat$fungene
  datako4<-subset(fungene,select=c(1:(nrow(group)+1)))

  datako4 %>%
    group_by(fungene) %>%
    summarise_all(sum) -> data

  data<-column_to_rownames(data,var = "fungene")
  data <- data*100
  data <- data %>% filter(apply(data,1,mean) > min.fun.abun)
  data <- t(data)
  data1 <- data.frame(data,group$group)
  colnames(data1) <- c(colnames(data),"Group")
  data1$Group <- as.factor(data1$Group)

  diff <- data1 %>%
    select_if(is.numeric) %>%
    map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')
  diff.dat<-diff
  if (padj) diff$p.value <- p.adjust(diff$p.value,"fdr")
  diff.dat$padj<-p.adjust(diff.dat$p.value,"fdr")
  if (filter) {
    diff1 <- diff %>% filter(p.value < 0.05)
    if (nrow(diff1) <= 0) {
      diff<-diff
    } else {
      diff<-diff1
    }
  }

  ## left barplot layout
  abun.bar <- data1[,c(diff$var,"Group")] %>%
    gather(variable,value,-Group) %>%
    group_by(variable,Group) %>%
    summarise(Mean = mean(value))

  kk1<-subset(fungene,select=c(fungene))
  kk1$variable<-kk1$fungene
  KKK<-merge(abun.bar,kk1,by="variable")
  abun.bar<-KKK
  abun.bar1<-abun.bar
  ## right point layout
  diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
  diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                              levels(data1$Group)[2]))
  diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]

  abun.bar$variable<-factor(abun.bar$variable,levels=rev(unique(abun.bar$variable)))

  diff.mean$fungene<-diff.mean$var
  ajj<-merge(diff.mean,kk1,by="fungene", all = TRUE,sort = FALSE)

  diff.mean1<-diff.mean

  diff.mean1[,"sig"]<-ifelse(diff.mean1$p.value<0.001,"***",
                             ifelse(diff.mean1$p.value<0.01,"**",
                                    ifelse(diff.mean1$p.value<0.05,"*","ns") ))

  cat("Totally seeking ",nrow(fungene)," functional genes.\n",sep = "")
  cat("Totally seeking ",nrow(diff.mean1)," functional genes, accounted for ",paste(sprintf("%0.2f",round(100*nrow(diff.mean1)/nrow(fungene))),"%",sep = ""),".\n",sep = "")

  file2=paste(export_path, "/tax4fun_FunGene","_",paste(tt[1],tt[2],sep = "-"),"_mean.txt",sep = "" )
  write.table(abun.bar,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Differential analysis results based on welch's t test using tax4fun prediction (",paste(tt[1],tt[2],sep = "-"),") has been exported to ","/",export_path,"",sep = "","\n")

  file2=paste(export_path, "/tax4fun_FunGene","_",paste(tt[1],tt[2],sep = "-"),"_onlydiffdata.txt",sep = "" )
  write.table(diff.mean1,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","Differential analysis results based on welch's t test using tax4fun prediction (",paste(tt[1],tt[2],sep = "-"),") has been exported to ","/",export_path,"",sep = "","\n")

  return(list(abun.bar=abun.bar,
              abun.bar1=abun.bar1,
              diff.mean=diff.mean,
              comparison=fp.dat$comparison_sel,
              otu_table=abun,
              diff.mean1=diff.mean1,
              diff.dat=diff.dat,
              sampledata=microchatTaxFun$qiimeotu$sampleda))
}



"plotDiffPathwayCombinedExtendbar" <- function(submchat,
                                               my_comparisons,
                                               xlabname,
                                               ko1.color=colorCustom(10,pal = "xiena"),
                                               add.bar.color=color_group,
                                               compare.color=color_group,
                                               ko1.select="Metabolism",
                                               padj=FALSE,
                                               min.fun.abun = 0.3,
                                               filter=TRUE,
                                               silva_path='/Volumes/BOOTCAMP/SILVA123',
                                               export_path) {

  xp<-list()
  for (i in 1:length(my_comparisons)) {
    microchatTaxFun<-calcTaxFun(submchat,
                                comparison_sel=paste(my_comparisons[[i]][1],
                                                     "-",
                                                     my_comparisons[[i]][2],sep = ""),
                                silva_path=silva_path,
                                export_path=export_path)

    microchatTaxFunDiff<-calcTaxFunDiff(ko1.select=ko1.select,     ###ko1.select="all",则使用所有ko1
                                        microchatTaxFun,
                                        padj=padj,
                                        min.fun.abun = min.fun.abun,   ###平均丰度大于该值才考虑组图
                                        filter=filter,          ###是否过滤掉差异不显著的
                                        export_path=export_path)

    ptreat2<-plotDiffPathwayExtendbar(microchatTaxFunDiff,
                                      ko2.bar.gradient=TRUE,
                                      gradient.num=20,
                                      xlabname=xlabname,
                                      pdf.width=10,
                                      ko2.text.size=4,
                                      layout.rt=c(0.4,5,4,4,1.5),
                                      ko1.color=ko1.color,
                                      add.bar.color=add.bar.color,
                                      compare.color=compare.color,
                                      siglabel="sig+label",
                                      export_path=export_path)

    xp[[i]]<-ptreat2
  }
  message("A list object has been exported. Patchwork package is recommended to be used for combined plots.")
  return(xp)
}



"plotDiffGeneVolcano" <- function(microchatTaxFun,
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
                                  export_path=export_path

) {

  fp.dat<-microchatTaxFun
  abun<-fp.dat$submchat$otu_table
  tt<-str_split_fixed(fp.dat$comparison_sel,'-',2)%>%as.character()
  control_treat<-strsplit(fp.dat$comparison_sel,"-")
  treat<-tt[2]
  control<-tt[1]

  select_comparison<-paste(tt[2],
                           tt[1],sep = "-")
  message("Control selected is ",tt[1])
  message("Treatment selected is ",tt[2])
  export_path<-paste(export_path,"/microbial functional prediction analysis/Functional gene/",tt[1],"-",tt[2],sep = "")
  dir.create(export_path, recursive = TRUE)
  group<-group_generate(abun)
  fungene<-fp.dat$fungene
  datako4<-fungene
  xlabname.sel<-xlabname[match(c(control,treat),unique(group$group))]

  datako4 %>%
    group_by(fungene) %>%
    summarise_all(sum) -> data

  data[,c("var","fungene")]<-str_split_fixed(data$fungene,";",2)
  data<-subset(data,select = -fungene)
  data<-column_to_rownames(data,var = "var")
  data <- data*1000000

  counts<-data%>%data.frame()
  group_order<-unique(substr(colnames(counts),start = 1,stop = 2))
  sampledata<-group%>%data.frame()
  sampledata<-sampledata[which(sampledata$group==tt[1]|sampledata$group==tt[2]),]
  my.counts.round<- counts %>% round()
  sampledata$group<-factor(sampledata$group,levels = unique(sampledata$group))
  group_list<-sampledata$group
  edger.counts <- edgeR::DGEList(counts = my.counts.round,
                                 group=factor(group_list),
                                 samples = sampledata)

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
  data<-result

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
      labs(title = paste(xlabname.sel[2]," regulation patterns compared to ",xlabname.sel[1],sep = ""))+
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
               family="serif",label = paste(xlabname.sel[2]," Down: ",label[1],sep = ""),
               min(DEG$logFC)*0.8, max((-log10(data$pvalue)))*0.95) +
      annotate('text',size = 4,color=sigcolor[3],
               family="serif", label = paste(xlabname.sel[2]," Up: ",label[3],sep = ""),
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
    sample.size=counts%>%colSums()
    if (length(unique(sample.size))==1) sample.size<-sample.size[1] else sample.size<-max(sample.size)
    otutab<-counts
    split_otux<-calc_indivGroup_abun(otutab)
    split_otux.mean<-(sapply(split_otux, rowMeans)/sample.size*100)%>%data.frame()
    split_otux.mean<-split_otux.mean%>%tibble::rownames_to_column(var = "name")

    abun.sel<-split_otux.mean[,c(1,match(control_treat[[1]],
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
      scale_x_discrete(labels=xlabname.sel)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      theme_void()+
      theme(text = element_text(family = "serif"),
            axis.text.x = element_text(colour='black',
                                       size=12,
                                       family = "serif",vjust = 1),
            legend.position = "none")


    barpv<-ggplot(node_phylum, aes(x = "", y = Freq, fill = Var1,alpha=Var1)) +
      geom_bar(stat = 'identity', width = 0.5) +
      coord_polar(theta = "y") +
      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.05,label=scales::percent(Freq/sum(node_phylum$Freq))),
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
      labs(title = paste(xlabname.sel[2]," regulation patterns compared to ",xlabname.sel[1],sep = ""))+
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
               family="serif",label = paste(xlabname.sel[2]," Down: ",labelx[1],sep = ""),
               -xcoord*0.9, max((-log10(data$pvalue)))*1.1) +
      annotate('text',size = 4,color=sigcolor[3],
               family="serif", label = paste(xlabname.sel[2]," Up: ",labelx[3],sep = ""),
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
    sample.size=counts%>%colSums()
    if (length(unique(sample.size))==1) sample.size<-sample.size[1] else sample.size<-max(sample.size)
    otutab<-counts
    split_otux<-calc_indivGroup_abun(otutab)
    split_otux.mean<-(sapply(split_otux, rowMeans)/sample.size*100)%>%data.frame()
    split_otux.mean<-split_otux.mean%>%tibble::rownames_to_column(var = "name")

    abun.sel<-split_otux.mean[,c(1,match(control_treat[[1]],
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
      scale_x_discrete(labels=xlabname.sel)+
      scale_alpha_manual(values = sigalpha,guide="none")+
      theme_void()+
      theme(text = element_text(family = "serif"),
            axis.text.x = element_text(colour='black',
                                       size=12,
                                       family = "serif",vjust = 1),
            legend.position = "none")


    barpv<-ggplot(node_phylum, aes(x = "", y = Freq, fill = Var1,alpha=Var1)) +
      geom_bar(stat = 'identity', width = 0.5) +
      coord_polar(theta = "y") +
      geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                    x= 1.05,label=scales::percent(Freq/sum(node_phylum$Freq))),
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

  DEG<-merge(DEG,abun.sel,by="name")
  #DEG$control<-control
  #DEG$treat<-treat
  file2=paste(export_path, "/single_volcano_",paste(control,treat,sep = "-"),"_oridata.txt",sep = "" )
  write.table(DEG,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  file2=paste(export_path, "/single_volcano_",paste(control,treat,sep = "-"),"_onlydiffdata.txt",sep = "" )
  write.table(DEG[which(DEG$change!="NOT"),],file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  cat("\n","Differential analysis results has been exported to“",export_path,"",sep = "")

  return(p2)
}



"plotDiffGeneCombinedVolcano" <- function(submchat,
                                          xlabname=xlabname,
                                          ncol=4,
                                          my_comparisons,
                                          remove.grid=TRUE,
                                          add.pieplot=TRUE,
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
                                          silva_path='/Volumes/BOOTCAMP/SILVA123',
                                          export_path) {
  xp<-list()
  i=1
  for (i in 1:length(my_comparisons)) {
    microchatTaxFun<-calcTaxFun(submchat,
                                comparison_sel=paste(my_comparisons[[i]][1],
                                                     "-",
                                                     my_comparisons[[i]][2],sep = ""),
                                silva_path=silva_path,
                                export_path=export_path)

    ### single volcano
    ptreat2<-plotDiffGeneVolcano(microchatTaxFun,
                                 xlabname=xlabname,
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
    xp[[i]]<-ptreat2
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
  ggsave(paste(export_path,"/microbial functional prediction analysis/Functional gene/Combined_volcano_",layout,".tiff",
               sep = ""),
         width = 6*width.sel,height = 6*height.sel,
         mp)
  ggsave(paste(export_path,"/microbial functional prediction analysis/Functional gene/Combined_volcano (RGB)_",layout,".pdf",
               sep = ""),colormodel="srgb",
         width = 6*width.sel,height = 6*height.sel,
         mp)
  ggsave(paste(export_path,"/microbial functional prediction analysis/Functional gene/Combined_volcano (CMYK)_",layout,".pdf",
               sep = ""),colormodel="cmyk",
         width = 6*width.sel,height = 6*height.sel,
         mp)

}
