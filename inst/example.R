####---- a simple workflow for microchat
#############-------1.data input、clean、filter & result output--------
#microchat::install.allpackage()
library(treeio)
library(microchat)
library(ape)
library(tidyverse)
library(ggforce)
library(igraph)
library(WGCNA)
x16s <- readRDS(system.file("extdata/x16s.rds", package = "microchat"))
abun<-x16s$otu_table
taxon<-x16s$taxon_table
tree11<-x16s$tree

abun<-otu_ultra_filter(abun, filter_num=3, seed = 1234)

####global setting
{
  export_path="hwj"                          ###relative path
  FileNamePrefix="param_"                 ###pre-suffix
  group.keep=NULL       
  xlabname=c("HP0","HP0.3","HP0.6","HP1.2")  
  my_comparisons <-list( c("CT","VH"),
                         c("CT","LH"),
                         c("CT","HH")
  )
  
  ###color
  color_taxa=adjustcolor(colorCustom(40,pal = "gygn"),alpha.f = 0.8)
  color_phylum=colorCustom(6,pal = "favor")  
  pal="favor"                                
  color_genus=colorCustom(6,pal = "favor")
  color_percent="black"
  color_group<-c(colorCustom(gnum,pal = "set1")[c(1,5)],colorCustom(gnum,pal = "set2")[1],colorCustom(gnum,pal = "set1")[3])
  color_group=adjustcolor(color_group,alpha.f = 0.6)
  color_NCM=colorCustom(3,pal = "favor")       
  color_BNTmodel=colorCustom(3,pal = "favor")  
  color_background = "grey90"
  color_circos=colorCustom(200,pal = "favor")
  sigcolor<-colorCustom(3,pal = "favor")
  pal_modulepreserve_taxa ="favor"
  pal_modulepreserve_cluster = "favor"
  link.pos="grey"       ###module preserve
  link.neg="#6E0000"    ###module preserve
  
  ###other
  linetype=1
  strictmod=TRUE
  method="anova"
  
  ###main theme
  main_theme = theme(panel.background=element_blank(),
                     panel.grid=element_blank(),
                     axis.line.x=element_line(size=0.5, colour="black"),
                     axis.line.y=element_line(size=0.5, colour="black"),
                     axis.ticks=element_line(color="black"),
                     axis.text=element_text(color="black", size=12),
                     legend.position="right",
                     legend.background=element_blank(),
                     legend.key=element_blank(),
                     legend.text= element_text(size=12),
                     text=element_text(family="serif", size=12),
                     aspect.ratio = 1)
  library(igraph)
  library(patchwork)
}

###mutiple phenotype files input
params<-microchatParamSink(FileNamePrefix="param_" )

mchat<-setMicrochat(abun,taxon,tree=tree11,params=params)

tidymchat<-tidyMicrochat(mchat,
                         group.keep=c("CT","VH","LH","HH"),  ### c("CT","SS","SA","SF","SC)
                         group.rm=NULL,    ### ct
                         sample.keep=NULL, ### c("ct1","ml1")
                         sample.rm=NULL,   ### c("ct1","ml1")
                         taxon.keep=NULL,  ### c("p__firm","c__bacter")
                         taxon.rm=NULL,    ### c("firm","cyan")
                         filter_pollution=TRUE,
                         filter_Archaea=TRUE)
submchat<-subsampleMicrochat(tidymchat,sample.size=NULL,export_path=export_path)

colSums(submchat$otu_table)%>%sum()
colSums(submchat$otu_table)
submchat$taxon_table$Phylum%>%unique()

#############-------2.microbial structure------

microchatcomobj<-calcMicrochatcom(submchat,
                                  taxa="Genus",
                                  taxa_num=10,
                                  export_path=export_path)

##bar
plotMicrocomBar(microchatcomobj,
                xlabname=xlabname,
                data_type="relative",
                sample_type="groups",
                chicklet=TRUE,
                color_taxa=color_taxa,
                color_group=color_group,
                barplot.barwidth=0.6,
                aspect.ratio=0.8,
                color_background = NA,
                color_border=NULL,
                export_path=export_path)

##pie
microchatcomobj<-calcMicrochatcom(submchat,
                                  taxa="Phylum",
                                  taxa_num=5,
                                  export_path=export_path)
microchatcomobj$plot_data$Taxo%>%unique()
plotMicrocomPie(microchatcomobj,
                xlabname=xlabname,
                data_type="relative",
                sample_type="groups",
                color_taxa=color_taxa,
                color_group=color_group,
                color_percent="black",
                export_path=export_path)

##heatmap
plotMicrocomHeatmap(microchatcomobj,
                    rescale=TRUE,
                    standardlization="0-1",
                    ###0-1 standardlization "scale", "center", "0-1"
                    data_type="relative",
                    sample_type="allsamples",
                    color_taxa=color_taxa,
                    color_group=color_group,
                    color.heatmap=c("white","black"),  
                    heatmap.add.line=TRUE,
                    heatmap.top.class=10,
                    heatmap.sample.angle=90,
                    heatmap.taxa.angle=0,
                    export_path=export_path)
##chord
par(family="serif")
plotMicrocomChord(microchatcomobj,
                  plot.reverse=FALSE,
                  text.size=1,
                  data_type="absolute",
                  sample_type="allsamples",
                  color_taxa=color_taxa,
                  color_group=color_group,
                  export_path=export_path)
##bubble
(plotMicrocomBubble(microchatcomobj,
                    data_type="absolute",
                    sample_type="allsamples",
                    color_taxa=color_taxa,
                    color_group=color_group,
                    bubble.color.taxa=TRUE,
                    bubble.max.size=9,
                    bubble.shape=19,    ## 15-20 solid, recommend; 0-14 hollow,
                    bubble.alpha=1,
                    aspect.ratio=0.4,
                    color_background = color_background,
                    color_border=NULL,
                    export_path=export_path)+
    theme(legend.position = "none"))/(plotMicrocomBubble(microchatcomobj,
                                                         data_type="absolute",
                                                         sample_type="groups",
                                                         color_taxa=color_taxa,
                                                         color_group=color_group,
                                                         bubble.color.taxa=TRUE,
                                                         bubble.max.size=9,
                                                         bubble.shape=19,    ## 15-20 solid, recommend; 0-14 hollow,
                                                         bubble.alpha=1,
                                                         aspect.ratio=0.4,
                                                         color_background = color_background,
                                                         color_border=NULL,
                                                         export_path=export_path)+
                                        theme(legend.position = "none"))

#############-------2.1 statistical analysis------
plotmicrochatComBoxplot(submchat,
                        taxa="Phylum",
                        taxa_num=10,
                        file.save=FALSE,
                        xlabname=xlabname,
                        ylim.fold=1.1, 
                        strictmod=TRUE,
                        method="t.test",
                        my_comparisons=my_comparisons,
                        color_group=color_group,
                        export_path=export_path) 

plotmicrochatComBarplot(submchat,
                        taxa="Phylum",
                        taxa_num=10,
                        file.save=FALSE,
                        xlabname=xlabname,
                        ylim.fold=1.1, 
                        strictmod=TRUE,
                        method="t.test",
                        my_comparisons=my_comparisons,
                        color_group=color_group,
                        export_path=export_path) 

#############-------2.2 differential analysi------
microchatTaxFunDiff<-calcmicrochatComDiff(submchat,
                                          comparison_sel="ct-xl",
                                          taxon.keep=NULL,    ###null 
                                          taxa_sel=c("Phylum","Class","Genus"),
                                          padj=FALSE,
                                          min.fun.abun = 1,   
                                          filter=FALSE)
plotTtestExtendBar(microchatTaxFunDiff,
                   abundance.order=FALSE,
                   abundance.order.rev=FALSE,
                   abundance.regulation.order=FALSE,
                   ko2.bar.gradient=TRUE,
                   ko2.bar.gradient.rev=FALSE,
                   ko1.color=colorCustom(10,pal = "xiena"),
                   add.bar.color=color_group,
                   compare.color=color_group,
                   gradient.num=20,
                   sig.text.size=3,
                   xlabname=xlabname,
                   sig.text.color="black",
                   ko2.text.size=4,
                   ko2.text.color="black",
                   add.bar.text.size=3,
                   add.bar.text.color="white",
                   layout.rt=c(1,4,4,4,1.5),
                   ko1.bar.alpha=1,
                   ko2.bar.alpha=1,
                   ko1.text.size=4,
                   ko1.text.color="black",
                   errorbar.shape=21, ### 22 23
                   diffbar.border.color = "black") 

par(mfrow=c(1,3))
groupname<-colnames(submchat$otu_table)%>%substr(start = 1,stop = 2)%>%unique()
for (gname in groupname) {
  plot_chord(submchat,
             all=FALSE,
             select_group=gname,###all=false
             select_taxa=c("Phylum","Genus"),
             transparency=0.4,
             directional=1,
             top_num=40,
             full_label=TRUE,
             char.start=1,
             char.stop=9,
             layout = "baseball",  #### regular  
             gridcol=color_taxa,
             diffHeight=0.06,
             annotationTrackHeight=0.03,
             link.arr.type=c("big.arrow"),
             small.gap = 1, big.gap = 5,
             textsize=0.4,
             lgdxpoi=1,lgdypoi=0,lgd_title="Taxa")

}

#############-------2.3.venn plot---------

microchatVennobj<-calcMicrochatVenn(submchat)

plotMicrochatAdvancedVenn(microchatVennobj,xlabname,
                          color_group ,
                          pieplot.mar = c(0.4,1.1,0.25,1),
                          linkplot.mar = c(-0.1,0.35,0.7,1.1),
                          export_path)

plotMicrochatVenn(microchatVennobj,
                  inter.point.color="black",
                  shade_color="grey75",
                  bar_color=color_group[3],
                  group_bar_color=color_group)

plotMicrochatVennDia(microchatVennobj)

plotMicrochatVennFlower(microchatVennobj,alpha_thres = 200) ###alpha_thres: max 255

###export to cytoscape
exportMicrochatVennCyto(microchatVennobj,
                        export_path=export_path)

#############-------3.assembly process----
###niche width
spec_gen<-calcMicrochatAssemblyNiche(submchat,
                                     export_path=export_path)
spec_gen$spe.gen.prop

plotMicrochatNicheWidth(spec_gen,
                        color_group=color_group,
                        export_path=export_path) 

plotMicrochatAssemblyNiche(spec_gen,nrow=1,
                           specified.thres=c(4,4,4,4,4),
                           sel.shape=c(21,19,23),
                           color_type=colorCustom(5,pal = "favor"),
                           export_path=export_path) 

###source of generalists and specialists
export_path0=paste(export_path,"/microbial assembly analysis/Generalist_Specialist Source",sep = "")
dir.create(export_path0, recursive = TRUE)
###layouts1
pdf(paste(export_path0,"/layout_Chord",".pdf",sep = ""),width = 30,height = 15)
plotMicrochatAssemblySourceChord(spec_gen,
                                 nrow=2,
                                 specified.thres=c(4,4,4,4,4),
                                 color_phylum=color_phylum,
                                 color_genus=color_genus)

dev.off()

###layouts2
require(igraph)
selected.layout<-c("circle","fr","nicely","kk","dh","tree","graphopt",
                   "gem","lgl","random","sphere","mds")
tt="circle"
for (tt in selected.layout) {
  dir.create(export_path0,recursive = TRUE)
  pdf(paste(export_path0,"/layout_",tt,".pdf",sep = ""),width = 30,height = 15)
  plotMicrochatAssemblySourceCircle(spec_gen,
                                    edge.curved=0.2,
                                    nrow=2,
                                    sel.layout=tt,  ### selected.layout[1]
                                    specified.thres=NULL,
                                    color_phylum=color_phylum,
                                    color_genus=color_genus)
  dev.off()
}



###neutral model
plot.ddat<-calcMicrochatAssemblyNCM(submchat,
                                    export_path=export_path)

dat<-sapply(plot.ddat, function(x){
  data.frame(x=x$Rsqr,
             y=1-x$Rsqr)})
rownames(dat)<-c("sto","deter")
dat

plotMicrochatAssemblyNCM(plot.ddat,
                         add.border=TRUE,  
                         color_background=NA, ## "grey95"  
                         color_pie.text="black",   
                         point.size=3,    
                         color_line=c("blue","blue","blue"), 
                         export_path=export_path)

#####BNT
library(picante)

BNT<-calcMicrochatAssemblyBNT(submchat,
                              bootstrap.thres=100,
                              nworker=4)

plotMicrochatAssemblyBNT(BNT,
                         geom="pie",  ## c("barplot","pie","chord")
                         pie.layout.nrow=1,
                         color_group=color_group,
                         color_model=colorCustom(5,pal = "favor"),
                         xlabname=xlabname,
                         export_path=export_path)->p1


#####NST--no tree input
library(ape)
library(iCAMP)
library(NST)

NST<-calcMicrochatAssemblyNST(submchat,
                              nworker=4,
                              dist.method="bray") 

plotMicrochatAssemblyNST(NST,
                         geom="pie",  ## c("barplot","pie","chord")
                         pie.layout.nrow=1,
                         color_group=color_group,
                         color_model=color_BNTmodel,
                         xlabname=xlabname,
                         export_path=export_path)->p2
p1/p2

#############-------4.diversity-------
###alpha-------
microchatAlphadivobj<-calcMicrochatAlphadiv(submchat,
                                            export_path=export_path)
microchatAlphadivobj$alphadiv
microchatAlphadivStatobj<-calcMicrochatAlphadivStat(microchatAlphadivobj,
                                                    alpha.index="Simpson",
                                                    strictmod=strictmod,
                                                    method=method,
                                                    comparison=my_comparisons,
                                                    export_path=export_path)
microchatAlphadivStatobj$alpha.stats
plotMicrochatAlphadiv(microchatAlphadivStatobj,
                      errorbar.line.add=FALSE,
                      errorbar.point.size=0.1,
                      y.point.adj=NULL,
                      seg=FALSE,
                      xlabname=xlabname,
                      color_group=color_group, 
                      color_background=color_background,
                      export_path=export_path)

####combinations
pp<-list()
alpha_index<-c("richness","Chao1" , "ACE", "Shannon" ,"Simpson" ,"pd","Pielou","Coverage")
for (t in alpha_index) {
  microchatAlphadivStatobj<-calcMicrochatAlphadivStat(microchatAlphadivobj,
                                                      alpha.index=t,
                                                      strictmod=strictmod,
                                                      method=method,
                                                      comparison=my_comparisons,
                                                      export_path=export_path)
  
  p<-plotMicrochatAlphadiv(microchatAlphadivStatobj,
                           errorbar.line.add=TRUE,
                           errorbar.point.size=0.1,
                           y.point.adj=NULL,
                           seg=TRUE,
                           xlabname=xlabname,
                           color_group=color_group,
                           color_background=NA,
                           export_path=export_path)+
    theme(panel.border = element_rect(fill = NA,color ="grey75",size = 1))
  
  pp[[t]]<-p
}


###beta-------
###permanova mrpp anoism
comm.pairwise.dat1<-pairwise.commdiff(submchat,
                                      plot.data=TRUE,  
                                      p.adjust.m = "BH",
                                      distance = 'jaccard',
                                      export_path=export_path)

comm.pairwise.dat2<-pairwise.commdiff(submchat,
                                      plot.data=TRUE, 
                                      p.adjust.m = "BH",
                                      distance = 'bray',
                                      export_path=export_path)

microchatBetadivobj1<-calcMicrochatBetadiv(submchat,
                                           distance = "bray", ###jaccard
                                           ordination="pcoa",
                                           export_path=export_path)

pb1<-plotMicrochatBetadiv(microchatBetadivobj1,
                          color_group=color_group,
                          color_background=NA,
                          xlabname=xlabname,
                          pcoa_3d=FALSE,
                          chull.shape=1,
                          text.size=2.5,
                          ellipse.linetype=2,
                          layout="ellipse",
                          export_path=export_path)

microchatBetadivobj2<-calcMicrochatBetadiv(submchat,
                                           distance = "jaccard",
                                           ordination="pcoa",
                                           export_path=export_path)

pb2<-plotMicrochatBetadiv(microchatBetadivobj2,
                          color_group=color_group,
                          color_background=NA,
                          xlabname=xlabname,
                          pcoa_3d=FALSE,
                          chull.shape=1,
                          text.size=2.5,
                          ellipse.linetype=2,
                          layout="ellipse",
                          export_path=export_path)

#############-------3. differential analysis-----------
#############-------5.1 lefse----------
##one-against-all
microchatLefseobj<-calcMicrochatLefse(submchat,
                                      lda_score=2,
                                      export_path=export_path)

lefsecladgram<-plotMicrochatLefse(microchatLefseobj,
                                  layout="radial",
                                  color_group=color_group,
                                  label_plot_display=5,  
                                  export_path=export_path)

##one-against-one
plotMicrochatComplexLefse(submchat,lda_score=2,
                          control="ct",
                          layout="radial",
                          label_plot_display=5,
                          color_group=color_group,
                          export_path=export_path)


plotMicrochatComplexLefse(submchat,lda_score=2,
                          control="CS",
                          layout="radial",
                          label_plot_display=5,
                          color_group=color_group,
                          export_path=export_path)
###
#############-------5.2edger&limma---------
plotMicrochatDiffCombinedVolcano(submchat,
                                 ncol=3,
                                 my_comparisons,
                                 xlabname=xlabname,
                                 remove.grid=TRUE,
                                 add.barplot=TRUE,
                                 sigcolor=colorCustom(5,pal = "favor"),
                                 sigalpha=c(1,0.4,1),
                                 layout="choice2",
                                 pvalue_thres=0.01,
                                 foldchange_thres=3,
                                 point.size=1,
                                 linetype="dashed",
                                 color_pie.text="black",
                                 color_background=NA,
                                 export_path) 

#############-------6.network analysis-----
####1.
rmnet_thres_calc<-rmnet_thres_calc(abun,
                                   filter_num = 2, 
                                   cor.method = "pearson") ####pearson spearman
###2.
rmnet_igraph<-rmnet_thres_select(trvalue,
                                 selected.thres=0.99,
                                 export_path=export_path)
rmnet_igraph$plot  ### preview
net<-rmnet_igraph$net  ### export matrix datasets for igraph object
rmt_mat<-lapply(net, function(x) {x$cleaned.matrix}) ### export matrix for igraph object
g<-rmnet_igraph$g   ### export igraph object
thres<-rmnet_igraph$result

###2.1 power-law fit
powerfit.r<-c()
for(i in 1:length(g)){
  g1 <- g[[i]]
  powerfit<-powerlaw_fit(g1)
  powerfit.r<-c(powerfit.r,powerfit$r2)
  cat("powerfit: ",powerfit$r2," p: ",powerfit$p,sep = "","\n")
}

powerlaw_logfit.r<-c()
for(i in 1:length(g)){
  g1 <- g[[i]]
  powerlaw_logfit<-powerlaw_logfit(g1)
  powerlaw_logfit.r<-c(powerlaw_logfit.r,powerlaw_logfit$r2)
  cat("powerlaw_logfit: ",powerlaw_logfit$r2," p: ",powerlaw_logfit$p,sep = "","\n")
}

###2.2 stability
###2.2.1 random removal
rmnet_cg<-rmnet_cytogephi_export(g,taxon,export_path =export_path)

randrm.data<-calcRMTnetRandrm(rmnet_igraph,rmnet_cg,split_otu,
                              #nperm=100,
                              export_path=export_path)

library(reshape2)
plotRMTnetRandrm(randrm.data,
                 geom="boxplot",  #barplot
                 select.rm.prop="50%",
                 strictmod=strictmod,
                 method=method,
                 comparison=my_comparisons,
                 xlabname=xlabname,
                 yaxis.italic=FALSE,
                 color_group=color_group,
                 export_path=export_path)

####2.2.2 taget removal
target.data<-calcRMTnetTargetrm(rmnet_igraph,rmnet_cg,split_otu,
                                nperm=100,
                                parallel=TRUE,
                                parallel.node.num=100,
                                #otu.sel=NULL,
                                export_path=export_path)

plotRMTnetTargrm(target.data,
                 geom="boxplot",
                 select.rm.modhub.num=7,
                 strictmod=strictmod,
                 method=method,
                 comparison=my_comparisons,
                 xlabname=xlabname,
                 yaxis.italic=FALSE,
                 color_group=color_group,
                 export_path=export_path)

####2.2.3 Vulner
Vulner<-calcRMTnetVulner(rmnet_igraph,rmnet_cg)

###3.adj matrix output
rmnet_matrix_export(net,export_path=export_path)

###4. graph preview
plot_rmnet(g,split_otu,thres,num=5,rows=1,layout = "fr",circle_plot = TRUE,
           vertex_size = 5,edge_size = 5,
           radius = 1.05,display_filter = FALSE)

###5.cytoscape&gephi file
rmnet_cg<-rmnet_cytogephi_export(g,taxon,export_path =export_path)

###5.1 module preserve
modpre_data<-calcRMTnetModPreserve(rmnet_cg,
                                   min.node.num=3,
                                   export_path=export_path)

##layout1
plotRMTnetModPreserveCircle(modpre_data,
                            edge.curved=1,
                            sel.layout="sphere",
                            ## "sphere"  "circle" "tree"
                            color_group=color_group,
                            color_cluster=colorCustom(13,pal = "xiena"))

##layout2
taxon<-submchat$taxon_table
g<-rmnet_igraph$g   ### export igraph object
plotRMTnetModPreserve(g,taxon,
                      pdf.sel=TRUE,
                      xlabname=xlabname,
                      bar.height=NULL, ### 
                      display.panel.border=FALSE,
                      display.interval.line=TRUE,
                      display.interval.line.style=c("grey","dotted",1),
                      pal_taxa =pal_modulepreserve_taxa,
                      pal_cluster = pal_modulepreserve_cluster,
                      color_bar.text="white",
                      min.node.num=6,
                      show.modulepairs.node=FALSE,
                      curvature=0.1,
                      link.pos=link.pos,  #"#3BBCC8"
                      link.neg=colorCustom(4,pal = pal_modulepreserve_taxa)[2],  #"#EC3822"
                      link.left.alpha=1,   
                      link.right.alpha=0.25,
                      link.pos.size=0.05,   
                      link.neg.size=0.2,   
                      link.clu.size=0.5,  
                      point.size=2,        
                      point.clu.size=20,  
                      export_path=export_path)

###5.2 positive & negative strength distribution
p1<-plot_RMTnetInter_Strength(rmnet_cg,
                              group.pos.coffe=0.1825,  
                              prop.coffe=1.05,         
                              bar.alpha=0.5,
                              xlabsname=xlabname,
                              add.line=FALSE,         
                              linetype=linetype,
                              bar.text.color="orange",
                              color_cor=c(link.pos,link.neg),
                              color_bar=c(link.pos,link.neg),
                              export_path=export_path)
###5.3 
p2<-plot_RMTnet_ModuleLink(rmnet_cg,
                           group.pos.coffe=0.1825,
                           prop.coffe=1.05,
                           bar.alpha=0.5,
                           xlabsname=xlabname,
                           add.line=FALSE,
                           linetype=linetype,
                           bar.text.color="orange",
                           color_cor=c(link.pos,link.neg),
                           color_bar=c(link.pos,link.neg),
                           export_path=export_path)

###5.4 node linked to the key nodes

sel.group=2 ### which group
nodes<-rmnet_cg$node_table
edges<-rmnet_cg$edge_table
node1s<-nodes[which(nodes$group==sel.group),]
node1s<-node1s[order(node1s$module,decreasing = FALSE),]
edge1s<-edges[which(edges$group==sel.group),]

plot_RMTnet_keynodeLink(node1s,edge1s,
                        keynode.link.display=FALSE, 
                        genus.coloring=TRUE,
                        same.point.size=FALSE, 
                        point.size=5,
                        curvature=0.2,
                        link.pos=link.pos,
                        link.neg=link.neg,
                        link.pos.size=1,
                        link.neg.size=1,
                        pal=pal,
                        text.size=2)

###5.5 node linked to the key nodes
plot_MutiRMTnet_keynodeLink(rmnet_cg,
                            nrows=1,
                            keynode.link.display=FALSE,
                            curvature=0.2,
                            link.pos="grey",
                            link.neg="#6E0000",
                            link.pos.size=0.3,
                            link.neg.size=0.3,
                            genus.coloring=TRUE,
                            pal="gygn",
                            same.point.size=FALSE,
                            point.size=1,
                            text.size=2,
                            export_path=export_path)


###6.perl-circos prefiles-----------
rmnet_circosfile(rmnet_cg,submchat$otu_table,export_path =export_path)

###7.perl-circos files generate---------
rmnet_circosToPerl(rmnet_cg,taxon,
                   export_path =export_path)

rmnet_circosColorToPerl(rmnet_cg,
                        color.circos=color_circos,
                        color.circos.alpha=NULL,
                        export_path =export_path)

###8.zipi plot、topological roles table-----------
rmnet_topo_compare(g,boot_times=999,
                   export_path =export_path)

plot_rmnet_ZPplot(rmnet_cg,abun,taxon,
                  type="mixed", # mixed all 都选
                  taxa="Phylum",
                  taxanum=10,
                  color_taxa=color_taxa,
                  export_path =export_path)

add.params=list(
  "Powerlaw R2 (log)"=powerfit.r
)

topo.dat<-rmnet_toporoles_export(rmnet_igraph,
                                 randcalc_times=100,  
                                 xlabname=xlabname,
                                 add.params=NULL,
                                 export_path =export_path)

###keynode table--------
keynode_tab<-rmnet_keynode_export(rmnet_cg,
                                  xlabname=xlabname,
                                  export_path =export_path)

###9. pei charts on cytoscape-style graphs-----------
plot_rmnet_cyto_pie_export(g,taxon,
                           color_text="#91091A",
                           color_taxa=color_taxa,
                           export_path =export_path)

plot_rmnet_modnode_barplot_export(g,taxon,
                                  geom.align="center", 
                                  point.direction="top",  #
                                  point.color.uniform=TRUE,  #
                                  point.density=20,  #
                                  point.panel.width=0.4,  ##
                                  point.size=0.4,  #
                                  point.text.size=10, #
                                  module.text.size=12, ##
                                  aspect.ratio=0.8, 
                                  color_module=color_taxa,
                                  color_group=color_group,
                                  export_path = export_path)

#############-------7.tree construction------
plotMicrochatTree(submchat,
                  open.angle=0,
                  equal.branch=FALSE,
                  layout = "rectangular", ###fan radial
                  # 'rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight' or 'ape'
                  pointsize=1,
                  branch.size=0.5,
                  color_taxa=c(color_taxa[1:5],rep("grey90",70)),
                  colors_group = color_group,
                  random = TRUE,geom.point = FALSE,
                  export_path=export_path)

plotMicrochatSingleTaxTree(submchat,
                           layout="fan",
                           color_taxa="darkgreen",
                           select.taxa="p__Firmicutes")

#############-------8. source tracking analysis------
microchatSourceobj<-calcMicrochatSource(submchat,
                                        nworker=4,
                                        sink_group = "CA",  
                                        source_group="CT") 

plot_sourcetrack(microchatSourceobj,
                 color_unknown="red",
                 sourcecolor = color_group,
                 plot_type="barplot",
                 repeatiton=FALSE,
                 barplot.nogreyback=FALSE,
                 export_path=export_path)

#############-------9. functional prediction------
####extended barplot
library(parallel)
microchatTaxFun<-calcTaxFun(submchat,
                            comparison_sel="ct-hl",
                            silva_path='/Volumes/BOOTCAMP/SILVA123',
                            export_path=export_path)

microchatTaxFunDiff<-calcTaxFunDiff(ko1.select="Metabolism",    
                                    microchatTaxFun,
                                    padj=FALSE,
                                    min.fun.abun = 0.3,   
                                    filter=TRUE,          
                                    export_path=export_path) 

plotDiffPathwayExtendbar(microchatTaxFunDiff,
                         ko2.bar.gradient=TRUE,
                         gradient.num=20,
                         xlabname=xlabname,
                         pdf.width=NULL,
                         ko2.text.size=4,
                         layout.rt=c(0.4,5,4,4,1.5),
                         ko1.color=colorCustom(10,pal = "gygn"),
                         add.bar.color=colorCustom(10,pal = "gygn"),
                         compare.color=color_group,
                         siglabel="sig+label",
                         export_path=export_path)

###functional genes
microchatTaxFunGeneDiff<-calcTaxFunGeneDiff(microchatTaxFun,
                                            padj=FALSE,
                                            min.fun.abun = 0.3,   
                                            filter=TRUE,        
                                            export_path)

plotDiffGeneExtendbar(microchatTaxFunGeneDiff,
                      xlabname=xlabname,
                      sig.text.size=3,
                      sig.text.color="black",
                      pdf.width=10,
                      layout.rt=c(4,4,1),
                      add.bar.color=colorCustom(50,pal = "gygn"),
                      compare.color=colorCustom(5,pal = "gygn"),
                      export_path=export_path)

plotDiffGeneVolcano(microchatTaxFun,
                    xlabname=xlabname,
                    remove.grid=TRUE,
                    add.pieplot=TRUE,
                    sigcolor=sigcolor,
                    sigalpha=c(1,0.4,1),
                    layout="choice2",
                    pvalue_thres=0.05,
                    foldchange_thres=1,
                    point.size=1,
                    linetype="dashed",
                    color_pie.text="white",
                    color_background=NA,
                    export_path=export_path)

plotDiffGeneCombinedVolcano(submchat,
                            xlabname=xlabname,
                            ncol=NULL,
                            my_comparisons,
                            remove.grid=TRUE,
                            add.pieplot=TRUE,
                            add.barplot=TRUE,
                            sigcolor=sigcolor,
                            sigalpha=c(1,0.4,1),
                            layout="choice2",
                            pvalue_thres=0.05,
                            foldchange_thres=1,
                            point.size=1,
                            linetype="dashed",
                            color_pie.text="black",
                            color_background=NA,
                            silva_path='/Volumes/BOOTCAMP/SILVA123',
                            export_path)

#############-------10.phenotype analysis--------
library(multcomp)
library(multcompView)
otu_table<-NULL
taxon_table<-NULL
tree<-NULL
export_path1<-paste(export_path,"/microbial parameteric analysis",sep = "")
color_background=NA
params<-microchatParamSink(FileNamePrefix="param_")
mchat<-setParamchat(otu_table,taxon_table,tree=tree,params=params)
tidymchat<-tidyMicrochat(mchat, group.keep=group.keep)
submchat<-subsampleMicrochat(tidymchat,sample.size=NULL,
                             export_path=export_path)

###10.1.
microchatParamobj<-calcMicrochatParam(submchat,
                                      paramfile.select="growth",
                                      export_path=export_path1)
##10.2
mcTable<-calcMicrochatParamTable(microchatParamobj,
                                 strictmod=TRUE,
                                 method="anova",
                                 comparison=my_comparisons,
                                 export_path=export_path)

##10.3
params.n<-calcMicrochat2Param(mcTable)
params.n1<-params.n[which(params.n$variable %in% mcTable$index),]
plotMicrochatParamBoxplot1(params.n1,color_group,label.x.angle=0)+
  ylab("Relative expression level")+
  theme(aspect.ratio = 0.8,legend.position = "none",
        axis.text.x = element_text(angle = 0,face = "bold.italic",hjust = 0.5),
        axis.line.x = element_line(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())

##10.4
### heatmap----
plotMicrochatParamHeatmap(microchatParamobj,
                          rescale=TRUE,
                          standardlization="normal",
                          ###0-1 standardlization "scale", "center", "0-1", "normal"
                          sample_type="groups",
                          color_taxa=color_taxa,
                          color_group=colorCustom(5,pal = "favor"),
                          color.heatmap=c("white",
                                          colorCustom(3,pal = "favor")[1]), 
                          heatmap.add.line=TRUE,
                          heatmap.top.class=3,
                          heatmap.sample.angle=0,
                          heatmap.taxa.angle=0,
                          export_path=export_path)
### boxplot----
plotMicrochatParamMutiBoxplot(microchatParamobj,
                              geom.line=FALSE,
                              xlabname=xlabname,
                              yaxis.italic=FALSE,
                              strictmod=FALSE,
                              panel.border=TRUE,
                              method="anova",
                              comparison=my_comparisons,
                              color_group=color_group,
                              color_backgroud=color_background,
                              export_path=export_path1)
### lineplot----
microchatParamStatobj<-calcMicrochatParamStat(microchatParamobj,
                                              select.index="MDA",
                                              strictmod=FALSE,
                                              method=method,
                                              comparison=my_comparisons,
                                              export_path=export_path1)
plotMicrochatParamLineplot(microchatParamStatobj,
                           xlabname=xlabname,
                           yaxis.italic=FALSE,
                           errorbar.line.add=TRUE,
                           errorbar.point.size=0.1,
                           y.point.adj=0.1,
                           seg=TRUE,
                           panel.border=TRUE,
                           add.spline=FALSE,
                           color_group=color_group,
                           color_backgroud=color_background,
                           export_path=export_path1)+ylim(0,1.2)

### barplot----
plotMicrochatParamMutiBarplot(microchatParamobj,
                              xlabname=xlabname,
                              yaxis.italic=TRUE,
                              strictmod=FALSE,
                              panel.border=TRUE,
                              method=method,
                              comparison=my_comparisons,
                              color_group=color_group,
                              color_backgroud=color_background,
                              export_path=export_path1)
### barplot----
plotMicrochatParamMuti2Barplot(microchatParamobj,
                               xlabname=xlabname,
                               yaxis.italic=FALSE,
                               panel.border=TRUE,
                               bar.border.size=1,
                               strictmod=TRUE,
                               method="anova",
                               comparison=my_comparisons,
                               color_group=color_group,
                               color_backgroud=color_background,
                               export_path=export_path1)

### bubble plot----
params<-read.delim('param_ss_liver_enzyme.txt',header = T,row.names = 1,stringsAsFactors = F)
class(tree11)
plotParametr(params[1:30,-c(1:3)],
             color_taxa=color_taxa,
             color_group=color_group,
             bubble.color.taxa=FALSE,
             bubble.max.size=9,
             bubble.shape=19,
             bubble.alpha=1,
             aspect.ratio=0.8,
             color_background=color_background,
             color_border=NULL,
             export_path=export_path1) 
params1<-read.delim('param_ss_liver_enzyme.txt',header = T,row.names = 1,stringsAsFactors = F)
plotParametr(params1,
             color_taxa=color_taxa,
             color_group=color_group,
             bubble.color.taxa=FALSE,
             bubble.max.size=9,
             bubble.shape=19,
             bubble.alpha=1,
             aspect.ratio=0.8,
             color_background=color_background,
             color_border=NULL,
             export_path=export_path1) 

#############-------10.5.data mining----
##net merge
g1<-mergeigraph(rmnet_igraph$g)
###obtain module eigengene
mes1<-get.moduleEigengenes(g1,submchat) 

###10.5.1
microchatEGHeatmap(mes1,axis.x.rev=TRUE,
                   color_sig=colorCustom(3,pal = "set2"),
                   method = "spearman") 

###10.5.2
mes<-mes1
microchatEGMatHeatmap(mes,NULL,axis.x.rev=TRUE,
                      color_sig=colorCustom(3,pal = "set2"),
                      method = "spearman") 
###10.5.3-------
library(randomForest)
rf.dat<-calcMicrochatRF(submchat)

opt.num<-get.RFThres(rf.dat,main_theme)

p.rf<-plotmicrochatRFbiomarker(opt.num,main_theme) 

( opt.num$p|p.rf$p.f)/
  ( opt.num$p1|p.rf$p.r)

otu.rf.f<-rf.dat$imp_otu$OTUid[1:opt.num$optimal.f]
otu.rf.r<-rf.dat$imp_otu$OTUid[1:opt.num$optimal.r]


###10.5.4-----
rmnet_igraph$plot
g1<-mergeigraph(rmnet_igraph$g)
rmnet_param_mult<-get.multparam.eigen(submchat,method = "spearman")
mes1<-get.moduleEigengenes(g1,submchat) 

mytheme=theme(panel.background = element_rect(fill = NA))
p1<-microchatEGParmHeatmap(rmnet_param_mult,
                           mes1,
                           add.anno=TRUE,
                           add.anno.range=1,
                           add.anno.alpha=0.8,
                           add.sec.param=TRUE,
                           mytheme=mytheme,
                           sig.text.size=3,
                           axis.x.rev=FALSE,
                           color_sig=colorCustom(3,pal = "set2"),
                           method = "spearman",
                           export_path)

mcbiomarker<-calcEGParmBiomarker(g1,submchat,
                                 rmnet_param_mult,
                                 cor.thres=0.2,
                                 p.thres=0.05,
                                 otu.mod.sel="ME3",
                                 param.mod.sel="gutjuction1")

rmnet_param_multxx<-get.param.eigen.t(submchat,rmnet_param_multx)
p3<-microchatEGParmHeatmap(rmnet_param_multxx,  ###或 rmnet_param_mult
                           ws$mes,
                           add.anno=TRUE,
                           add.anno.range=0,
                           add.anno.alpha=0.8,
                           add.sec.param=FALSE,
                           mytheme=mytheme,
                           sig.text.size=3,
                           axis.x.rev=FALSE,
                           color_sig=sigcolor,
                           method = "spearman",
                           export_path)
