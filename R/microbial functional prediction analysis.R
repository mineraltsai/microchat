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
  library(Tax4Fun)
  data(kegg_anno)
  kegg_anno=kegg_anno
  dir.create(export_path, recursive = TRUE)

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
  ggsave(paste(export_path,"/",paste(id[1],id[2],sep = "-"),".pdf",sep = ""),pp)
  cat("图片已保存至相对路径","/",export_path,"","下",sep = "")
  return(list(extendbar=pp,pred_res=funpred))
}



"readtaxdata" <- function(submchat,comparison = "M0-M5")
{

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

  select_comparison<-strsplit(comparison,"-")
  select_comparison<-select_comparison[[1]]

  slee<-match(sampleda$group,select_comparison)%>%as.numeric()
  pos<-which(!is.na(slee))

  abun_sel<-abun[,pos]
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

"tax4funpred" <- function(QIIMESingleData,
                          kegg_anno,
                          folderReferenceData = '/Volumes/BOOTCAMP/SILVA123',
                          export_path='~/Desktop/文章撰写/jbh1/tax4fun') {

  #KEGG 功能基因（KO 第四级）丰度预测
  Tax4FunOutput <- Tax4Fun::Tax4Fun(QIIMESingleData, folderReferenceData, fctProfiling = TRUE, refProfile = 'UProC', shortReadMode = TRUE, normCopyNo = TRUE)
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
  Tax4FunOutput <- Tax4Fun::Tax4Fun(QIIMESingleData, folderReferenceData, fctProfiling = FALSE, refProfile = 'UProC', shortReadMode = TRUE, normCopyNo = TRUE)
  tax4fun_pathway <- as.data.frame(t(Tax4FunOutput$Tax4FunProfile))
  pathway <- rownames(tax4fun_pathway)
  write.table(cbind(pathway, tax4fun_pathway),
              paste(export_path,"/pathway.",groupname.use,".txt",sep = ""),
              row.names = FALSE, sep = '\t', quote = FALSE)

  ##为 KEGG 代谢通路映射 KO 分类
  tax4fun_pathway$KO3_id <- t(data.frame(strsplit(rownames(tax4fun_pathway), ';')))[ ,1]
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
    scale_x_discrete(limits = levels(diff.mean$var)) +
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
