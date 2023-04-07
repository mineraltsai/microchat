"microchatalpha" <- function(otu, tree = NULL, base = 2) {
  otu<-t(otu)
  est <- vegan::estimateR(otu)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- vegan::diversity(otu, index = 'shannon', base = base)
  Simpson <- vegan::diversity(otu, index = 'simpson')    #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(otu == 1) / rowSums(otu)

  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- picante::pd(otu, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  return(result)
}

"calcMicrochatAlphadiv" <- function(submchat,
                                    export_path="microbial diversity analysis") {

  dir.create(export_path, recursive = TRUE)

  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  otu<-submchat$otu_table
  tree<-submchat$tree
  alphadiv<-microchatalpha(otu, tree, base = 2)

  colnames(alphadiv)<-c("Richness", "Shannon", "Simpson", "Pielou", "Chao1", "ACE", "Coverage", "Pd" )
  alpha.table<-tibble::rownames_to_column(alphadiv,var = "sample")
  alpha.table$group<-substr(alpha.table$sample,start = 1,stop = 2)
  alpha.table<-subset(alpha.table, select = c(which(colnames(alpha.table)=="group"),
                                              which(colnames(alpha.table)!="group")))

  write.table(alpha.table,file = paste(export_path,"/alpha diversity.txt",sep = ""),row.names = FALSE,quote = FALSE, sep = "\t")
  cat("--------------------------------------------------------------------------------------")

  alpha.table$group<-factor(alpha.table$group,levels = unique(alpha.table$group))

  ###normality test
  nor.data<-subset(alpha.table,select = -c(group,sample))
  nptest<-data.frame()
  for (i in 1:length(colnames(nor.data))) {
    sw<-norm.test(nor.data[,i])
    {
      nptest<-rbind(nptest,sw)
    }
  }

  cat("--------------------------------------------------------------------------------------")
  ###Homogeneity of variance
  variance.data<-subset(alpha.table,select = -c(sample))
  vttest<-data.frame()
  for (i in 2:length(colnames(variance.data))) {
    input.data=subset(variance.data,select = c(group,i))
    index=colnames(input.data)[2]
    vt<-variance.test(input.data,index,alpha = 0.05)
    {
      vttest<-rbind(vttest,vt)
    }
  }


  vttest[,c("index","none","char")]<-str_split_fixed(vttest$data.name," ",3)
  nptest$index<-vttest$index

  alphadiv<-alpha.table
  alphadiv$group<-factor(alphadiv$group,levels = unique(alphadiv$group))

  ###output statistics table
  all.index<-colnames(alphadiv)[3:length(colnames(alphadiv))]
  cat("\n Please use one of",all.index," in the plot function !!!\n")

  summ_fcbv.all<-data.frame()
  for (index in all.index) {
    # calculate the mean and standard error of each group according to the alpha index
    summ_fcbv<- alphadiv %>%
      group_by(group) %>%
      plyr::summarise(
        mean = aggregate(as.formula(paste(index, " ~ group",sep = "")),data=alphadiv,FUN = "mean"),
        sd = aggregate(as.formula(paste(index, " ~ group",sep = "")),data=alphadiv,FUN = "sd")
      ) %>%data.frame()

    summ_fcbv1<-cbind(summ_fcbv$mean,summ_fcbv$sd)
    summ_fcbv1<-subset(summ_fcbv1, select=c(-3))
    colnames(summ_fcbv1)<-c("group","mean","sd")

    summ_fcbv1<-summ_fcbv1[1:length(unique(alphadiv$group)),]
    summ_fcbv1$group<-factor(summ_fcbv1$group,levels = unique(alphadiv$group))

    gpnum<-data.frame()
    for (gname in unique(alphadiv$group)) {
      gnum<-which(alphadiv$group==gname)%>%length()
      newd<-data.frame(group=gname,
                 n=gnum,
                 index=index)
      {
        gpnum<-rbind(gpnum,newd)
      }
    }

    summ_fcbv1<-left_join(summ_fcbv1,gpnum,by="group")
    summ_fcbv1$se<-summ_fcbv1$sd/sqrt(summ_fcbv1$n)
    summ_fcbv1$group<-factor(summ_fcbv1$group,levels = unique(summ_fcbv1$group))

    {
      summ_fcbv.all<-rbind(summ_fcbv.all,summ_fcbv1)
    }
  }

  write.table(summ_fcbv.all,file = paste(export_path,"/alpha diversity.stat.txt",sep = ""),row.names = FALSE,quote = FALSE, sep = "\t")
  message("\n","The alpha diversity was used for statistical analysis. You could check it.","\n")

  message("\n","The alpha diversity has been statistically analyzed. You could check it.","\n")

  microchatAlphadivobj<-list(alpha.table,nptest,vttest,summ_fcbv.all,all.index,all.group)

  names(microchatAlphadivobj)[1]<-"alphadiv"
  names(microchatAlphadivobj)[2]<-"np.test"
  names(microchatAlphadivobj)[3]<-"vt.test"
  names(microchatAlphadivobj)[4]<-"statistics"
  names(microchatAlphadivobj)[5]<-"all.index"
  names(microchatAlphadivobj)[6]<-"all.group"

  class(microchatAlphadivobj) <- c("microchat","data.frame")

  return(microchatAlphadivobj)
}


"calcMicrochatAlphadivStat" <- function(microchatAlphadivobj,
                                        alpha.index="shannon",
                                        strictmod=FALSE,
                                        method=c("t.test","wilcox.test","anova","kruskal.test"),
                                        comparison=my_comparisons,
                                        export_path="microbial diversity analysis") {


  method<-match.arg(method)
  dir.create(export_path, recursive = TRUE)

  if (class(microchatAlphadivobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  all.group<-microchatAlphadivobj$all.group

  alphadiv<-microchatAlphadivobj$alphadiv
  alphadiv$group<-factor(alphadiv$group,levels = unique(alphadiv$group))
  ###output alpha table
  alphadiv.use<-reshape::melt(alphadiv)

  ### words match
  all.index<-microchatAlphadivobj$all.index
  index<-all.index[match(tolower(alpha.index),tolower(all.index))]

  ###define new alpha table according to the selected alpha index
  data_poi<-alphadiv.use[which(alphadiv.use$variable==index),]
  data_poi$group <- factor(data_poi$group , levels = unique(data_poi$group))

  ###define new statistics table according to the selected alpha index
  data_err<-microchatAlphadivobj$statistics
  data_err<-data_err[which(data_err$index==index),]
  data_err$group <- factor(data_err$group , levels = unique(data_err$group))


  ###assume the sample number varied
  groupnum<-unique(alphadiv$group)%>%length()
  samplenum<-unique(alphadiv$sample)%>%length()

  #my_comparisons<-calcComparisonOrder(my_comparisons,data_poi)

  if (!strictmod) {

    ##!strictmod
    ###assume the variance test of all data
    np.test<-microchatAlphadivobj$np.test
    vp.test<-microchatAlphadivobj$vt.test

    matchnp.par<-np.test[np.test$p.value>0.05,]%>%rownames()%>%as.numeric()
    matchvp.par<-vp.test[vp.test$p.value>0.05,]%>%rownames()%>%as.numeric()

    match.num<-union(matchnp.par,matchvp.par)[order(union(matchnp.par,matchvp.par),
                                                    decreasing = FALSE)]

    diff.num<-setdiff(1:length(np.test$index),match.num)%>%length()
    diff.order<-setdiff(1:length(np.test$index),match.num)%>%as.numeric()
    diff.order<-diff.order[order(diff.order,decreasing = FALSE)]
    diff.index<-all.index[diff.order]
    match.index<-setdiff(all.index,diff.index)

    if (method %in% c("t.test","wilcox.test")) {

      if (index %in% match.index) {
        method= "t.test"
      } else {
        method= "wilcox.test"
      }

      if ( method== "t.test") {
        sig_label_new=NULL
       sta<-lapply(my_comparisons, function(x){
        alphadiv.use.selected1<-data_poi[
          which(data_poi$group==x[1] | data_poi$group==x[2]),]
        if (length(which(data_poi$group==x[1]))==length(which(data_poi$group==x[2]))) {
          fit<-t.test(value~group, alphadiv.use.selected1,
                      paired = TRUE, alternative = 'two.sided')
        } else {
          fit<-t.test(value~group, alphadiv.use.selected1,
                      paired = FALSE, alternative = 'two.sided')
        }
        stats<-data.frame(
          index=index,
          group1=x[1],
          group2=x[2],
          mean1=mean(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[1])]),
          mean2=mean(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[2])]),
          max1=max(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[1])]),
          max2=max(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[2])]),
          statistics=fit$statistic,
          p.value=fit$p.value,
          method=fit$method)
      })

      data.alpha<-data.frame()
      for (tt in 1:length(sta)) {
        data.al<-sta[[tt]]
        {
          data.alpha<-rbind(data.alpha,data.al)
        }
      }


      data.alpha$sig<-ifelse(data.alpha$p.value<0.001,"***",
                             ifelse(data.alpha$p.value<0.01,"**",
                                    ifelse(data.alpha$p.value<0.05,"*","ns")))

      data.alpha$sd1<-data_err$sd[match(data.alpha$group1, data_err$group)]
      data.alpha$sd2<-data_err$sd[match(data.alpha$group2, data_err$group)]
      data.alpha$yvalue1<-data.alpha$mean1+data.alpha$sd1
      data.alpha$yvalue2<-data.alpha$mean2+data.alpha$sd2

      {group.selected<-c(unique(data.alpha$group1),unique(data.alpha$group2))

        value.selected<-c(unique(data.alpha$yvalue1),unique(data.alpha$yvalue2))

        max.value.selected<-value.selected%>%max()

        max.value.group.selected<-group.selected[match(max.value.selected,
                                                       value.selected)]


      }

      y.ratio<-(max(data_poi$value)-min(data_poi$value))/30
      data.alpha$order<-(match(data.alpha$group2,
                               data_err$group)-match(data.alpha$group1,
                                                     data_err$group))

      data.alpha$error<-(data.alpha$max1-data.alpha$max2)%>%abs()
      data.alpha<-data.alpha[order(data.alpha$order,data.alpha$error,decreasing = FALSE),]
}

      if ( method== "wilcox.test") {
        sig_label_new=NULL
        sta<-lapply(my_comparisons, function(x){
        alphadiv.use.selected1<-data_poi[
          which(data_poi$group==x[1] | data_poi$group==x[2]),]
        if (length(which(data_poi$group==x[1]))==length(which(data_poi$group==x[2]))) {
          fit<-wilcox.test(value~group, alphadiv.use.selected1,
                      paired = TRUE, alternative = 'two.sided')
        } else {
          fit<-wilcox.test(value~group, alphadiv.use.selected1,
                      paired = FALSE, alternative = 'two.sided')
        }
        stats<-data.frame(
          index=index,
          group1=x[1],
          group2=x[2],
          mean1=mean(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[1])]),
          mean2=mean(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[2])]),
          max1=max(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[1])]),
          max2=max(alphadiv.use.selected1$value[
            which(alphadiv.use.selected1$group==x[2])]),
          statistics=fit$statistic,
          p.value=fit$p.value,
          method=fit$method)
      })

      data.alpha<-data.frame()
      for (tt in 1:length(sta)) {
        data.al<-sta[[tt]]
        {
          data.alpha<-rbind(data.alpha,data.al)
        }
      }


      data.alpha$sig<-ifelse(data.alpha$p.value<0.001,"***",
                             ifelse(data.alpha$p.value<0.01,"**",
                                    ifelse(data.alpha$p.value<0.05,"*","ns")))

      data.alpha$sd1<-data_err$sd[match(data.alpha$group1, data_err$group)]
      data.alpha$sd2<-data_err$sd[match(data.alpha$group2, data_err$group)]
      data.alpha$yvalue1<-data.alpha$mean1+data.alpha$sd1
      data.alpha$yvalue2<-data.alpha$mean2+data.alpha$sd2

      {group.selected<-c(unique(data.alpha$group1),unique(data.alpha$group2))

        value.selected<-c(unique(data.alpha$yvalue1),unique(data.alpha$yvalue2))

        max.value.selected<-value.selected%>%max()

        max.value.group.selected<-group.selected[match(max.value.selected,
                                                       value.selected)]


      }

      y.ratio<-(max(data_poi$value)-min(data_poi$value))/30
      data.alpha$order<-(match(data.alpha$group2,
                               data_err$group)-match(data.alpha$group1,
                                                     data_err$group))

      data.alpha$error<-(data.alpha$max1-data.alpha$max2)%>%abs()
      data.alpha<-data.alpha[order(data.alpha$order,data.alpha$error,decreasing = FALSE),]
      }

      #data.alpha$error<-(data.alpha$max1-data.alpha$max2)%>%abs()
      #data.alpha<-data.alpha[order(data.alpha$order),]
    }

    if (method %in% c("anova","kruskal.test")) {
      variance.data<-subset(alphadiv,select=c(-sample))
      ###anova
      aovtest<-data.frame()
      test.b <- c()
      test.c<- data.frame()
      if (index %in% match.index) {
        method= "anova"
      } else {
        method= "kruskal.test"
      }

      if ( method == "anova") {
        variance.data<-subset(variance.data, select=c(1,match(match.index,colnames(variance.data))))
        for (i in 2:length(colnames(variance.data))) {
          if (i==2) {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-aov.test(input.data,colnames(input.data)[2])
            tuk1<-multcomp::glht(fit1,linfct=multcomp::mcp(group="Tukey"))
            res1 <- multcomp::cld(tuk1,alpah=0.05)
            ano<-summary(fit1)%>%'[['(1)%>%as.data.frame()
            ano1<-ano$`Pr(>F)`%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            aovtest<-rbind(aovtest,ano1)
            test.b <- cbind(test.b,res1$mcletters$Letters)
          } else {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-aov.test(input.data,colnames(input.data)[2])
            tuk1<-multcomp::glht(fit1,linfct=multcomp::mcp(group="Tukey"))
            res1 <- multcomp::cld(tuk1,alpah=0.05)
            ano<-summary(fit1)%>%'[['(1)%>%as.data.frame()
            ano1<-ano$`Pr(>F)`%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            {
              aovtest<-cbind(aovtest,ano1)
              test.b <- cbind(test.b,res1$mcletters$Letters)
            }
          }
        }


        test.b<-test.b%>%data.frame()
        colnames(test.b)<-colnames(variance.data)[2:ncol(variance.data)]

        sig_label<-test.b
        sig_label$group<-rownames(sig_label)
        index_order<-which(colnames(sig_label)==index)
        sig_label<-subset(sig_label,select = c(index_order,group))
        sig_label$group <- factor(sig_label$group , levels = unique(sig_label$group))

        data_max<-subset(data_poi,select=c(1,4))%>%
          group_by(group)%>%
          summarise_all(max)

        sig_label_new<-merge(sig_label,data_max,by="group")
        colnames(sig_label_new)[2]<-"alpha"
        sig_label_new$group <- factor(sig_label_new$group , levels = unique(sig_label$group))
        sig_label_new<-merge(sig_label_new,data_err,by="group")
        sig_label_new$valuey<-sig_label_new$value+sig_label_new$se

        x = c(1:length(data_poi$sample))
        y = data_poi$value

        ##方差分析
        fit1 <- aov(value ~  group,
                    data = data_poi)
        ##事后比较
        tuk1<-glht(fit1,linfct=mcp(group="Tukey"))
        ##显著性标志
        res1 <- cld(tuk1,alpah=0.05)
        ##显著性数值
        ano<-summary(fit1)%>%'[['(1)%>%as.data.frame()
        pvalue<-ano$`Pr(>F)`[1]
        if (pvalue<0.001) pvalue=0.001

        data.alpha<-data.frame(
          index=index,
          p.value=pvalue,
          method=method
        )
        y.ratio=NULL
        }

      if ( method == "kruskal.test") {
        variance.data<-subset(variance.data, select=c(1,match(diff.index,colnames(variance.data))))

        input.data<-reshape2::melt(variance.data)
        input.data1<-kruskalmuticomp.t(input.data)
        for (i in 2:length(colnames(variance.data))) {
          if (i==2) {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-kruskal(input.data,colnames(input.data)[2])
            tuk1<-dunn.test(input.data,colnames(input.data)[2])
            ano<-fit1$p.value
            ano1<-ano%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            aovtest<-rbind(aovtest,ano1)
            test.c<-rbind(test.c,tuk1)
          } else {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-kruskal(input.data,colnames(input.data)[2])
            tuk1<-dunn.test(input.data,colnames(input.data)[2])
            ano<-fit1$p.value
            ano1<-ano%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            {
              aovtest<-cbind(aovtest,ano1)
              test.c<-rbind(test.c,tuk1)
            }
          }
        }

        input.data<-input.data1[which(input.data1$variable==index),]
        sig_label_new<-input.data
        samplenum<-variance.data$group%>%table()%>%max()
        sig_label_new$valuey<-sig_label_new$mean+sig_label_new$std*sqrt(samplenum)
        sig_label_new$alpha<-sig_label_new$Letters


        aovtest=aovtest[which(colnames(aovtest)==index)]
        if (aovtest<0.001) aovtest=0.001
        p.value<-aovtest
        data.alpha<-data.frame(
          index=index,
          p.value=p.value,
          method=method
        )
        y.ratio=NULL
      }

    }

  } else {

    ##strictmod
    ###assume the variance test of all data
    np.test<-microchatAlphadivobj$np.test
    vp.test<-microchatAlphadivobj$vt.test

    matchnp.par<-np.test[np.test$p.value>0.05,]%>%rownames()%>%as.numeric()
    matchvp.par<-vp.test[vp.test$p.value>0.05,]%>%rownames()%>%as.numeric()

    match.num<-intersect(matchnp.par,matchvp.par)[order(intersect(matchnp.par,matchvp.par),
                                                        decreasing = FALSE)]

    diff.num<-setdiff(1:length(np.test$index),match.num)

    diff.order<-setdiff(1:length(np.test$index),match.num)%>%as.numeric()
    diff.order<-diff.order[order(diff.order,decreasing = FALSE)]
    diff.index<-all.index[diff.order]
    match.index<-setdiff(all.index,diff.index)

    if (method %in% c("t.test","wilcox.test")) {

      if (index %in% match.index) {
        method= "t.test"
      } else {
        method= "wilcox.test"
      }

      if ( method== "t.test") {
        sig_label_new=NULL
        sta<-lapply(my_comparisons, function(x){
          alphadiv.use.selected1<-data_poi[
            which(data_poi$group==x[1] | data_poi$group==x[2]),]
          if (length(which(data_poi$group==x[1]))==length(which(data_poi$group==x[2]))) {
            fit<-t.test(value~group, alphadiv.use.selected1,
                        paired = TRUE, alternative = 'two.sided')
          } else {
            fit<-t.test(value~group, alphadiv.use.selected1,
                        paired = FALSE, alternative = 'two.sided')
          }
          stats<-data.frame(
            index=index,
            group1=x[1],
            group2=x[2],
            mean1=mean(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[1])]),
            mean2=mean(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[2])]),
            max1=max(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[1])]),
            max2=max(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[2])]),
            statistics=fit$statistic,
            p.value=fit$p.value,
            method=fit$method)
        })

        data.alpha<-data.frame()
        for (tt in 1:length(sta)) {
          data.al<-sta[[tt]]
          {
            data.alpha<-rbind(data.alpha,data.al)
          }
        }


        data.alpha$sig<-ifelse(data.alpha$p.value<0.001,"***",
                               ifelse(data.alpha$p.value<0.01,"**",
                                      ifelse(data.alpha$p.value<0.05,"*","ns")))

        data.alpha$sd1<-data_err$sd[match(data.alpha$group1, data_err$group)]
        data.alpha$sd2<-data_err$sd[match(data.alpha$group2, data_err$group)]
        data.alpha$yvalue1<-data.alpha$mean1+data.alpha$sd1
        data.alpha$yvalue2<-data.alpha$mean2+data.alpha$sd2

        {group.selected<-c(unique(data.alpha$group1),unique(data.alpha$group2))

          value.selected<-c(unique(data.alpha$yvalue1),unique(data.alpha$yvalue2))

          max.value.selected<-value.selected%>%max()

          max.value.group.selected<-group.selected[match(max.value.selected,
                                                         value.selected)]


        }

        y.ratio<-(max(data_poi$value)-min(data_poi$value))/30
        data.alpha$order<-(match(data.alpha$group2,
                                 data_err$group)-match(data.alpha$group1,
                                                       data_err$group))

        data.alpha$error<-(data.alpha$max1-data.alpha$max2)%>%abs()
        data.alpha<-data.alpha[order(data.alpha$order,data.alpha$error,decreasing = FALSE),]
      }

      if ( method== "wilcox.test") {
        sig_label_new=NULL
        sta<-lapply(my_comparisons, function(x){
          alphadiv.use.selected1<-data_poi[
            which(data_poi$group==x[1] | data_poi$group==x[2]),]
          if (length(which(data_poi$group==x[1]))==length(which(data_poi$group==x[2]))) {
            fit<-wilcox.test(value~group, alphadiv.use.selected1,
                             paired = TRUE, alternative = 'two.sided')
          } else {
            fit<-wilcox.test(value~group, alphadiv.use.selected1,
                             paired = FALSE, alternative = 'two.sided')
          }
          stats<-data.frame(
            index=index,
            group1=x[1],
            group2=x[2],
            mean1=mean(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[1])]),
            mean2=mean(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[2])]),
            max1=max(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[1])]),
            max2=max(alphadiv.use.selected1$value[
              which(alphadiv.use.selected1$group==x[2])]),
            statistics=fit$statistic,
            p.value=fit$p.value,
            method=fit$method)
        })

        data.alpha<-data.frame()
        for (tt in 1:length(sta)) {
          data.al<-sta[[tt]]
          {
            data.alpha<-rbind(data.alpha,data.al)
          }
        }


        data.alpha$sig<-ifelse(data.alpha$p.value<0.001,"***",
                               ifelse(data.alpha$p.value<0.01,"**",
                                      ifelse(data.alpha$p.value<0.05,"*","ns")))

        data.alpha$sd1<-data_err$sd[match(data.alpha$group1, data_err$group)]
        data.alpha$sd2<-data_err$sd[match(data.alpha$group2, data_err$group)]
        data.alpha$yvalue1<-data.alpha$mean1+data.alpha$sd1
        data.alpha$yvalue2<-data.alpha$mean2+data.alpha$sd2

        {group.selected<-c(unique(data.alpha$group1),unique(data.alpha$group2))

          value.selected<-c(unique(data.alpha$yvalue1),unique(data.alpha$yvalue2))

          max.value.selected<-value.selected%>%max()

          max.value.group.selected<-group.selected[match(max.value.selected,
                                                         value.selected)]


        }

        y.ratio<-(max(data_poi$value)-min(data_poi$value))/30
        data.alpha$order<-(match(data.alpha$group2,
                                 data_err$group)-match(data.alpha$group1,
                                                       data_err$group))

        data.alpha$error<-(data.alpha$max1-data.alpha$max2)%>%abs()
        data.alpha<-data.alpha[order(data.alpha$order,data.alpha$error,decreasing = FALSE),]
      }

      #data.alpha$error<-(data.alpha$max1-data.alpha$max2)%>%abs()
      #data.alpha<-data.alpha[order(data.alpha$order),]
    }

    if (method %in% c("anova","kruskal.test")) {
      variance.data<-subset(alphadiv,select=c(-sample))
      ###anova
      aovtest<-data.frame()
      test.b <- c()
      test.c<- data.frame()
      if (index %in% match.index) {
        method= "anova"
      } else {
        method= "kruskal.test"
      }

      if ( method == "anova") {
        variance.data<-subset(variance.data, select=c(1,match(match.index,colnames(variance.data))))
        for (i in 2:length(colnames(variance.data))) {
          if (i==2) {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-aov.test(input.data,colnames(input.data)[2])
            tuk1<-multcomp::glht(fit1,linfct=multcomp::mcp(group="Tukey"))
            res1 <- multcomp::cld(tuk1,alpah=0.05)
            ano<-summary(fit1)%>%'[['(1)%>%as.data.frame()
            ano1<-ano$`Pr(>F)`%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            aovtest<-rbind(aovtest,ano1)
            test.b <- cbind(test.b,res1$mcletters$Letters)
          } else {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-aov.test(input.data,colnames(input.data)[2])
            tuk1<-multcomp::glht(fit1,linfct=multcomp::mcp(group="Tukey"))
            res1 <- multcomp::cld(tuk1,alpah=0.05)
            ano<-summary(fit1)%>%'[['(1)%>%as.data.frame()
            ano1<-ano$`Pr(>F)`%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            {
              aovtest<-cbind(aovtest,ano1)
              test.b <- cbind(test.b,res1$mcletters$Letters)
            }
          }
        }


        test.b<-test.b%>%data.frame()
        colnames(test.b)<-colnames(variance.data)[2:ncol(variance.data)]

        sig_label<-test.b
        sig_label$group<-rownames(sig_label)
        index_order<-which(colnames(sig_label)==index)
        sig_label<-subset(sig_label,select = c(index_order,group))
        sig_label$group <- factor(sig_label$group , levels = unique(sig_label$group))

        data_max<-subset(data_poi,select=c(1,4))%>%
          group_by(group)%>%
          summarise_all(max)

        sig_label_new<-merge(sig_label,data_max,by="group")
        colnames(sig_label_new)[2]<-"alpha"
        sig_label_new$group <- factor(sig_label_new$group , levels = unique(sig_label$group))
        sig_label_new<-merge(sig_label_new,data_err,by="group")
        sig_label_new$valuey<-sig_label_new$value+sig_label_new$se

        x = c(1:length(data_poi$sample))
        y = data_poi$value

        ##方差分析
        fit1 <- aov(value ~  group,
                    data = data_poi)
        ##事后比较
        tuk1<-glht(fit1,linfct=mcp(group="Tukey"))
        ##显著性标志
        res1 <- cld(tuk1,alpah=0.05)
        ##显著性数值
        ano<-summary(fit1)%>%'[['(1)%>%as.data.frame()
        pvalue<-ano$`Pr(>F)`[1]
        if (pvalue<0.001) pvalue=0.001

        data.alpha<-data.frame(
          index=index,
          p.value=pvalue,
          method=method
        )
        y.ratio=NULL
      }

      if ( method == "kruskal.test") {
        variance.data<-subset(variance.data, select=c(1,match(diff.index,colnames(variance.data))))

        input.data<-reshape2::melt(variance.data)
        input.data1<-kruskalmuticomp.t(input.data)
        for (i in 2:length(colnames(variance.data))) {
          if (i==2) {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-kruskal(input.data,colnames(input.data)[2])
            tuk1<-dunn.test(input.data,colnames(input.data)[2])
            ano<-fit1$p.value
            ano1<-ano%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            aovtest<-rbind(aovtest,ano1)
            test.c<-rbind(test.c,tuk1)
          } else {
            input.data=subset(variance.data,select = c(group,i))
            fit1<-kruskal(input.data,colnames(input.data)[2])
            tuk1<-dunn.test(input.data,colnames(input.data)[2])
            ano<-fit1$p.value
            ano1<-ano%>%data.frame()
            colnames(ano1)<-colnames(variance.data)[i]
            {
              aovtest<-cbind(aovtest,ano1)
              test.c<-rbind(test.c,tuk1)
            }
          }
        }

        input.data<-input.data1[which(input.data1$variable==index),]
        sig_label_new<-input.data
        samplenum<-variance.data$group%>%table()%>%max()
        sig_label_new$valuey<-sig_label_new$mean+sig_label_new$std*sqrt(samplenum)
        sig_label_new$alpha<-sig_label_new$Letters


        aovtest=aovtest[which(colnames(aovtest)==index)]
        if (aovtest<0.001) aovtest=0.001
        p.value<-aovtest
        data.alpha<-data.frame(
          index=index,
          p.value=p.value,
          method=method
        )
        y.ratio=NULL
      }

    }

  }

  #write.table(data.alpha,file = paste(export_path,"/alpha diversity.stat.diff.txt",sep = ""),row.names = FALSE,quote = FALSE, sep = "\t")
  message("\n","The alpha diversity has been statistically analyzed. You could check it.","\n")

  microchatAlphadivStatobj<-list(data.alpha,data_poi,data_err,index,method,y.ratio,sig_label_new,all.group)

  names(microchatAlphadivStatobj)[1]<-"alpha.stats"
  names(microchatAlphadivStatobj)[2]<-"data_poi"
  names(microchatAlphadivStatobj)[3]<-"data_err"
  names(microchatAlphadivStatobj)[4]<-"index"
  names(microchatAlphadivStatobj)[5]<-"method"
  names(microchatAlphadivStatobj)[6]<-"y.ratio"
  names(microchatAlphadivStatobj)[7]<-"sig"
  names(microchatAlphadivStatobj)[8]<-"all.group"

  class(microchatAlphadivStatobj) <- c("microchat","data.frame")

  return(microchatAlphadivStatobj)
}


"plotMicrochatAlphadiv" <- function(microchatAlphadivStatobj,
                                    errorbar.pos.adj=FALSE,
                                    errorbar.line.add=FALSE,
                                    seg=TRUE,
                                    errorbar.point.size=0.1,
                                    y.point.adj=NULL,
                                    color_group=colorCustom(5,pal = "gygn"),
                                    color_backgroud="grey90",
                                    export_path="microbial diversity analysis") {

  dir.create(export_path, recursive = TRUE)

  if (class(microchatAlphadivStatobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data_poi<-microchatAlphadivStatobj$data_poi
  data_err<-microchatAlphadivStatobj$data_err
  data.alpha<-microchatAlphadivStatobj$alpha.stats
  index<-microchatAlphadivStatobj$index
  if (is.null(y.point.adj)) {
    y.ratio<-microchatAlphadivStatobj$y.ratio
  } else {
    y.ratio<-microchatAlphadivStatobj$y.ratio*(1+y.point.adj)
    }
  method<-microchatAlphadivStatobj$method
  sig_label_new<-microchatAlphadivStatobj$sig


  all.group<-microchatAlphadivStatobj$all.group
  sel.group<-unique(data_poi$group)

  sel.order<-match(sel.group,all.group)

  colors<-color_group[sel.order]
  names(colors)<-unique(data_poi$group)


 if (method %in% c("t.test","wilcox.test")) {
   if (errorbar.pos.adj)
    for (kk in 1:nrow(data.alpha)) {
      if (data.alpha[kk,]$yvalue1>data.alpha[kk,]$max1){
        data.alpha[kk,]$yvalue1<-data.alpha[kk,]$yvalue1

      } else {
        data.alpha[kk,]$yvalue1<-data.alpha[kk,]$max1
      }


      if (data.alpha[kk,]$yvalue2>data.alpha[kk,]$max2){
        data.alpha[kk,]$yvalue2<-data.alpha[kk,]$yvalue2

      } else {
        data.alpha[kk,]$yvalue2<-data.alpha[kk,]$max2
      }
    }

  p<-ggplot(data=data_poi,
            aes(x = group,y = value,fill = group)) +
    geom_errorbar(data = data_err,
                  aes(x = group, y = mean, group = index, ymin = mean-sd, ymax = mean+sd),
                  colour = "grey50",width=.25)+
    geom_boxplot(outlier.shape = NA,
                 width = 0.5,
                 color = "white")+
    geom_jitter(data=data_poi,size=3,alpha=0.5,
                aes(group=index,color=group))+
    geom_point(aes(group=index,color=group))


  if (!seg) data.seg<-calcSegmentOrder(data.alpha,y.ratio,data_poi)
  if (seg) data.seg<-calcSegment2Order(data.alpha,y.ratio,data_poi)

  for (tt in 1:nrow(data.seg)) {
    x1<-data.seg$order1[tt]
    x2<-data.seg$order2[tt]

    y1.use<-data.seg$y1[tt]
    y2.use<-data.seg$y2[tt]
    y3.use<-data.seg$y3[tt]
    y4.use<-data.seg$y4[tt]

    if (y2.use>y3.use) {
      y3.use<-y2.use+y.ratio
    }

    if (y1.use>y4.use) {
      y4.use<-y1.use+y.ratio
    }

    if (y3.use>y4.use) {
      y3.use<-y3.use
    } else {
      y3.use<-y4.use
    }

    label.use=data.seg$sig[tt]

    p<-p+
      annotate("pointrange", color="grey50",size=errorbar.point.size,
               x = x1, y = y1.use, ymin = y1.use, ymax = y1.use)+
      annotate("pointrange", color="grey50",size=errorbar.point.size,
               x = x2, y = y2.use, ymin = y2.use, ymax = y2.use)+

      annotate("segment",color="grey50", x = x1, y = y1.use, xend = x1, yend = y3.use)+
      annotate("segment",color="grey50", x = x2, y = y2.use, xend = x2, yend = y3.use)+
      annotate("segment",color="grey50", x = x1, y = y3.use, xend = x2, yend = y3.use)+
      annotate("text",family="serif",x = (x1+x2)/2, y = y3.use+y.ratio/2,label=label.use)

    if (errorbar.line.add) p<-p+annotate("segment", color="grey50",x = x1-0.05, y = y1.use, xend = x1+0.05, yend = y1.use)+
      annotate("segment",color="grey50", x = x2-0.05, y = y2.use, xend = x2+0.05, yend = y2.use)

  }


  p<-p+ylab(index) +
    theme(#panel.border = element_blank(),
      axis.ticks.length = unit(0.2,"lines"),
      axis.ticks = element_line(color='black'),
      panel.background = element_rect(fill = color_backgroud),
      #axis.line = element_line(colour = "black"),
      axis.title.x=element_blank(),
      axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
      axis.text.y=element_text(colour='black',size=10,family = "serif"),
      axis.text.x=element_text(colour = "black",size = 10,
                               angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
      #strip.background =  element_blank(),
      legend.position = "none",aspect.ratio = 1)


  p<-p+scale_discrete_manual(values=colors,
                             aesthetics = "colour")+
    scale_discrete_manual(values=colors,
                          aesthetics = "fill")

  p<-p+annotate("text",family="serif",
                x = 1, y = min(data_poi$value),label=method)

 }

 if (method %in% c("anova","kruskal.test")){

    if (method=="anova") {p<-ggplot() +
      geom_errorbar(data = data_err,
                    aes(x = group, y = mean, group = index,
                        ymin = mean-sd, ymax = mean+sd),
                    colour = "grey50",width=.25)+
      geom_boxplot(data=data_poi,
                   aes(x = group,y = value,fill = group),
                   outlier.shape = NA,
                   width = 0.5,
                   color = "white")+
      geom_jitter(data=data_poi,size=3,alpha=0.5,
                  aes(x=group, y=value, group=index,color=group))+
      geom_point(data=data_poi,
                 aes(x=group, y=value, group=index,color=group)) +
      geom_text(data = sig_label_new,
                aes(x = group,y = valuey,label = alpha),
                size = 5,color = "black",family = "serif") +
      #theme_test()+
      ylab(index) +
      theme(#panel.border = element_blank(),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_line(color='black'),
        panel.background = element_rect(fill = color_backgroud),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=10,family = "serif"),
        axis.text.x=element_text(colour = "black",size = 10,
                                 angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
        legend.position = "none",aspect.ratio = 1)

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")
    p<-p +
      annotate('text', label = sprintf('italic(pvalue) == %.3f', signif(data.alpha$p.value,3)),
               x = 1, y = max(data_poi$value)*1.02,
               size = 3, parse = TRUE,family = "serif")

    p<-p+annotate("text",family="serif",
                  x = 1, y = min(data_poi$value),label=method)}

   if (method=="kruskal.test") {
     p.value<-sig_label_new

  p<-ggplot() +
     geom_errorbar(data = data_err,
                   aes(x = group, y = mean, group = index,
                       ymin = mean-sd, ymax = mean+sd),
                   colour = "grey50",width=.25)+
     geom_boxplot(data=data_poi,
                  aes(x = group,y = value,fill = group),
                  outlier.shape = NA,
                  width = 0.5,
                  color = "white")+
     geom_jitter(data=data_poi,size=3,alpha=0.5,
                 aes(x=group, y=value, group=index,color=group))+
     geom_point(data=data_poi,
                aes(x=group, y=value, group=index,color=group)) +
     #theme_test()+
     ylab(index) +
     theme(#panel.border = element_blank(),
       axis.ticks.length = unit(0.2,"lines"),
       axis.ticks = element_line(color='black'),
       panel.background = element_rect(fill = color_backgroud),
       #axis.line = element_line(colour = "black"),
       axis.title.x=element_blank(),
       axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
       axis.text.y=element_text(colour='black',size=10,family = "serif"),
       axis.text.x=element_text(colour = "black",size = 10,
                                angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
       legend.position = "none",aspect.ratio = 1)

   p<-p+scale_discrete_manual(values=colors,
                              aesthetics = "colour")+
     scale_discrete_manual(values=colors,
                           aesthetics = "fill")
   p<-p +
     annotate('text', label = sprintf('italic(pvalue) == %.3f', signif(p.value,3)),
              x = 1, y = max(data_poi$value)*1.02,
              size = 3, parse = TRUE,family = "serif")

   p<-p+annotate("text",family="serif",
                 x = 1, y = min(data_poi$value),label=method)}

  }

  ggsave(paste(export_path,"/alpha_diversity (",index,") .pdf",sep = ""),p)
  cat("\n","Alpha diversity barplot has been exported. Please check it.","\n")

  return(p)

}


"calcMicrochatBetadiv" <- function(submchat,
                                   distance = "bray",
                                   ordination="PCoA",
                                   export_path="microbial diversity analysis") {


  dir.create(export_path, recursive = TRUE)

  message(" bray distance based on abudance table，jaccard distance based on community structure.")
  message(" distance: bray, jaccard, wei_unifrac or unwei_unifrac")

  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  tax <- submchat$taxon_table
  otu <- submchat$otu_table
  sel.group<-unique(substr(colnames(otu),start = 1,stop = 2))

  otu <- data.frame(t(otu))
  abun<-submchat$otu_table
  group<-data.frame(group=substr(colnames(abun),start = 1, stop = 2),
                    sampleid=colnames(abun))
  rownames(group)<-group$sampleid
  group$group<-factor(group$group,levels = unique(group$group))

  adonis_drug <- vegan::adonis2(otu~group, data=group,
                         distance = distance, permutations = 999)

  bray_dis <- vegdist(otu, method = distance)
  if (ordination=="pcoa") pcoa <- cmdscale(bray_dis, k = (nrow(otu) - 1), eig = TRUE)
  if (ordination=="nmds") pcoa<-metaMDS(bray_dis,k = 2)
  pcoa_exp <- pcoa$eig/sum(pcoa$eig)
  pcoa_exp<-round(pcoa_exp*100,2)

  ####pvalue & r2
  pvalue<-adonis_drug$`Pr(>F)`[1]
  r2<-adonis_drug$R2[1]

  ####ordination
  scores<-pcoa$points[,1:3]%>%data.frame()
  colnames(scores)<-paste("PCoA",1:ncol(scores),sep = "")
  eig<-pcoa_exp[1:3]


  scores$site<-substr(rownames(scores),start = 1, stop = 2)
  scores$site<-factor(scores$site,levels = unique(scores$site))

  site<-cbind(group,scores)
  site$group<-factor(site$group,levels = unique(site$group))

  tax <- submchat$taxon_table
  otu <- submchat$otu_table
  otu <- data.frame(t(otu))
  group<-data.frame(group=substr(colnames(abun),start = 1, stop = 2),
                    sampleid=colnames(abun))
  rownames(group)<-group$sampleid
  group$group<-factor(group$group,levels = unique(group$group))

  library(vegan)
  adonis_drug <- vegan::adonis2(otu~group, data=group,
                         distance = distance, permutations = 999)

  bray_dis <- vegdist(otu, method = distance)
  if (ordination=="pcoa") pcoa <- cmdscale(bray_dis, k = (nrow(otu) - 1), eig = TRUE)
  if (ordination=="nmds") pcoa<-metaMDS(bray_dis,k = 2)
  pcoa_exp <- pcoa$eig/sum(pcoa$eig)
  pcoa_exp<-round(pcoa_exp*100,2)

  dune.pairwise.adonis <- pairwiseAdonis::pairwise.adonis(x=otu, factors=group$group, sim.function = "vegdist",
                                          sim.method = distance,
                                          p.adjust.m = "BH",
                                          reduce = NULL,
                                          perm = 999)

  file2=paste(export_path,"/beta_div (",ordination,"-",distance,")_beta_scores.txt",sep = "")
  write.table(scores,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity based on Permanova (",ordination,"-",distance,") ",",has been exported to","/",export_path,"",sep = "","\n")

  file2=paste(export_path,"/beta_div (",ordination,"-",distance,")_beta_exp.txt",sep = "")
  write.table(eig,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity based on Permanova_proportion ",",has been exported to","/",export_path,"",sep = "","\n")

  file2=paste(export_path,"/beta_div(",ordination,"-",distance,")_eig.txt",sep = "")
  write.table(pcoa_exp,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity based on Permanova_eig ",",has been exported to ","/",export_path,"",sep = "","\n")

  file2=paste(export_path,"/beta_div(",ordination,"-",distance,")_permanova_stat.txt",sep = "")
  write.table(dune.pairwise.adonis,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity results based on Permanova ",",has been exported to ","/",export_path,"",sep = "","\n")


  microchatBetadivobj<-list(scores,
                            eig,
                            pcoa_exp,
                            dune.pairwise.adonis,
                            site,
                            pvalue,
                            r2,
                            ordination,
                            distance,
                            all.group,
                            sel.group)

  names(microchatBetadivobj)[1]<-"scores"
  names(microchatBetadivobj)[2]<-"eig"
  names(microchatBetadivobj)[3]<-"pcoa_exp"
  names(microchatBetadivobj)[4]<-"dune.pairwise.adonis"
  names(microchatBetadivobj)[5]<-"site"
  names(microchatBetadivobj)[6]<-"pvalue"
  names(microchatBetadivobj)[7]<-"r2"
  names(microchatBetadivobj)[8]<-"ordination"
  names(microchatBetadivobj)[9]<-"distance"
  names(microchatBetadivobj)[10]<-"all.group"
  names(microchatBetadivobj)[11]<-"sel.group"

  class(microchatBetadivobj) <- c("microchat","data.frame")
  return(microchatBetadivobj)
}


"plotMicrochatBetadiv" <- function(microchatBetadivobj,
                                   color_group=colorCustom(5,pal = "gygn"),
                                   color_backgroud="grey90",
                                   pcoa_3d=FALSE,
                                   chull.shape=1,
                                   text.size=2.5,
                                   ellipse.linetype=2,
                                   layout=c("ellipse","chull","ellipse+chull"),
                                   export_path="microbial diversity analysis") {

  layout<-match.arg(layout)
  dir.create(export_path, recursive = TRUE)

  if (class(microchatBetadivobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  site<-microchatBetadivobj$site
  scores<-microchatBetadivobj$scores
  eig<-microchatBetadivobj$eig
  pcoa_exp<-microchatBetadivobj$pcoa_exp
  dune.pairwise.adonis<-microchatBetadivobj$dune.pairwise.adonis
  pvalue<-microchatBetadivobj$pvalue
  r2<-microchatBetadivobj$r2

  timer<-data.frame()
  for (tt in unique(scores$site)) {
    time.a <- scores[scores$site == tt,
    ][chull(scores[scores$site == tt,
                   c("PCoA1", "PCoA2")]),]
    {
      timer<-rbind(timer,time.a)
    }

  }

  if (length(color_group)<length(unique(site$group))) {
    stop("please provide more colors for coloring groups !!!")
  } else {
    all.group<-microchatBetadivobj$all.group
    sel.group<-microchatBetadivobj$sel.group
    sel.order<-match(sel.group,all.group)
    color.use<-color_group[sel.order]
    names(color.use)<-unique(site$group)
  }

  p <- ggplot(data = site,aes(x = PCoA1, y = PCoA2, color = group)) +
    geom_point(aes(x = PCoA1, y = PCoA2, color = group), size = 2)



  if (layout== "ellipse")
    p<-p+stat_ellipse(aes(x = PCoA1, y = PCoA2, color = group), level = 0.95,
                      linetype = ellipse.linetype, show.legend = FALSE)

  if (layout== "chull")
    p<-p+ggalt::geom_encircle(aes(fill=site), alpha = 0.1, show.legend = F,
                              expand=0,spread=0,s_shape=chull.shape)+
    geom_polygon(data=timer,aes(x=PCoA1,y=PCoA2,color=site,fill=site),
                 linetype=ellipse.linetype,alpha=0.1,show.legend = F)

  if (layout=="ellipse+chull")
    p<-p+stat_ellipse(aes(x = PCoA1, y = PCoA2, color = group), level = 0.95,
                      linetype = ellipse.linetype, show.legend = FALSE)+
    ggalt::geom_encircle(aes(fill=site), alpha = 0.1, show.legend = F,
                         expand=0,spread=0.5,s_shape=chull.shape) +
    geom_polygon(data=timer,aes(x=PCoA1,y=PCoA2,color=site,fill=site),
                 linetype=ellipse.linetype,alpha=0.1,show.legend = F)



  p<-p+
    ggrepel::geom_text_repel(aes(x = PCoA1, y = PCoA2, label = sampleid),

                             size = text.size,family="serif",
                             box.padding = unit(0.3, 'lines'), show.legend = FALSE,
                             max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
    scale_color_manual(values = color.use) +
    labs(x = "pcoa1", y = "pcoa2", color = '')


  p<-p+

    scale_discrete_manual(values=color.use,
                          aesthetics = "fill")+
    xlab(paste("PCoA1: ",pcoa_exp[1],"%",sep = ""))+
    ylab(paste("PCoA2: ",pcoa_exp[2],"%",sep = ""))+
    theme(panel.background = element_rect(fill = color_backgroud),
          text = element_text(family = "serif"),aspect.ratio = 1)



  colors<-color.use
  colors_3d<-colors
  names(colors)<-unique(scores$site)
  ypos<-scores$PCoA2%>%max()
  xpos<-scores$PCoA1%>%min()


  ypostop<-scores$PCoA2%>%max()
  yposbot<-scores$PCoA2%>%min()

  xposrig<-scores$PCoA1%>%max()
  xposlef<-scores$PCoA1%>%min()

  yposi<-ifelse(abs(ypostop)>abs(yposbot),abs(ypostop),abs(yposbot))
  xposi<-ifelse(abs(xposrig)>abs(xposlef),abs(xposrig),abs(xposlef))

  xposi<-(xposi*1.05*10)%>%as.integer()
  xposi<-xposi+1
  xposi<-xposi/10

  yposi<-(yposi*1.05*10)%>%as.integer()
  yposi<-yposi+1
  yposi<-yposi/10

  p<-p+
    annotate('text', label = sprintf('italic(pvalue) == %.3f', signif(pvalue,3)),
             x = xpos*0.9, y =ypos*1.05, size = 3, parse = TRUE,family = "serif")+
    annotate('text', label = sprintf('italic(r2) == %.3f', signif(r2,3)),
             x = xpos*0.5, y =ypos*1.05, size = 3, parse = TRUE,family = "serif")+
    annotate('text', label = "permanova",
             x = xpos*0.9, y =ypos*1.1, size = 5, parse = TRUE,family = "serif")

  if (layout== "chull") p<-p+xlim(-xposi,xposi)+ylim(-yposi,yposi)

  if (pcoa_3d) {
    scores_3d<-scores
    for (k in 1:3) {
      colnames(scores_3d)[k]<-paste(colnames(scores)[k],": ",eig[k],"%",sep = "")
    }
    colors_3d <- colors_3d[as.numeric(as.factor(site$group))]

    s3d<-scatterplot3d::scatterplot3d(scores_3d[,1:3],
                                      pch = 16,       # 点形状
                                      color=colors_3d,   # 点颜色
                                      cex.symbols = 2     )

    legend("top",
           legend = unique(site$group),
           col =  colors,
           pch = 16,
           inset = -0.1,
           xpd = TRUE,
           horiz = TRUE)

    text(s3d$xyz.convert(scores_3d[,c(1,2,3)] + 2),
         labels = rownames(scores_3d),
         cex = 0.8,col = "black")


  }
  ordination<-microchatBetadivobj$ordination
  distance<-microchatBetadivobj$distance

  ggsave(paste(export_path,"/beta_diversity (",ordination,") based on",distance," .pdf",sep = ""),p)
  cat("\n","Beta diversity cluster plot has been exported. Please check it.","\n")
  return(p)
}



"dune.pairwise.anosim"<- function(otu,
                                  p.adjust.m= "BH",
                                  distance = 'bray',
                                  export_path="cs2/microbial diversity analysis"){

  dir.create(paste(export_path,"/beta diversity/anosim",sep = ""), recursive = TRUE)
  group<-group_generate(otu)
  otu <- data.frame(t(otu), stringsAsFactors = FALSE)
  otu <- otu[rowSums(otu[])>0,]

  sample_id <- row.names(otu)
  part_group  <- data.frame(group,stringsAsFactors = FALSE)
  #print(group)

  group <- part_group[which (part_group[,1] %in% sample_id),]
  anosim_result_otu <- vegan::anosim(otu, group$group, permutations = 999,distance =distance)

  group_name <- unique(group$group)
  anosim_result_two <- NULL
  for (i in 1:(length(group_name) - 1)) {
    for (j in (i + 1):length(group_name)) {
      group_ij <- subset(group, group %in% c(group_name[i], group_name[j]))
      otu_ij <- otu[group_ij$sample, ]
      anosim_result_otu_ij <- anosim(otu_ij,
                                     group_ij$group,
                                     permutations = 999,
                                     distance = distance)

      if (anosim_result_otu_ij$signif <= 0.001) Sig <- '***'
      else if (anosim_result_otu_ij$signif <= 0.01) Sig <- '**'
      else if (anosim_result_otu_ij$signif <= 0.05) Sig <- '*'
      else Sig <- ""
      anosim_result_otu_ij$statistic<-sprintf('%.3f',round(anosim_result_otu_ij$statistic,3))
      anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep = ' vs '), distance,
                                                      anosim_result_otu_ij$statistic, anosim_result_otu_ij$signif, Sig))

    }

  }
  anosim_result_two<-anosim_result_two%>%data.frame()
  colnames(anosim_result_two)<-c("comparison","distance","statistic (r)","pvalue","sig")
  anosim_result_two$p.adj =p.adjust(anosim_result_two$pvalue,method=p.adjust.m)

  for (tt in 1:nrow(anosim_result_two)) {
    if (anosim_result_two$p.adj[tt] <= 0.001) anosim_result_two$sig[tt] <- '***'
    else if (anosim_result_two$p.adj[tt] <= 0.01) anosim_result_two$sig[tt] <- '**'
    else if (anosim_result_two$p.adj[tt] <= 0.05) anosim_result_two$sig[tt] <- '*'
    else anosim_result_two$sig[tt] <- ""
  }
  anosim_result_two<-subset(anosim_result_two,select=c(1,2,3,4,6,5))

  if (distance == 'bray') distance<-"Bray-curtis"
  if (distance == 'jaccard') distance<-"Jaccard"

  anosim_result_two$distance<-distance
  anosim_result_two$p.adj<-sprintf('%.3f',round(anosim_result_two$p.adj,3))

  file2=paste(export_path,"/beta diversity/anosim/anosim_",distance,".txt",sep = "")
  write.table(anosim_result_two,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity based on ANOSIM_",distance," has been exported to","/",export_path,"",sep = "","\n")

  anosim_result_two<-subset(anosim_result_two,select=c(1,2,3,5,6))
  colnames(anosim_result_two)[4]<-"pvalue"
  return(anosim_result_two)
}

"dune.pairwise.adonis"<- function(otu,
                                  p.adjust.m = "BH",
                                  distance = 'bray',
                                  export_path="cs2/microbial diversity analysis") {
  dir.create(paste(export_path,"/beta diversity/adonis",sep = ""), recursive = TRUE)
  group<-group_generate(otu)
  otux <- data.frame(t(otu), stringsAsFactors = FALSE)
  otux <- otux[rowSums(otux[])>0,]

  group$group<-factor(group$group,levels = unique(group$group))


  dune.pairwise.adonis <- pairwiseAdonis::pairwise.adonis(x=otux,
                                          factors=group$group,
                                          sim.function = "vegdist",
                                          sim.method = distance,
                                          p.adjust.m = p.adjust.m,
                                          reduce = NULL,
                                          perm = 999)
  dune.pairwise.adonis$F.Model<-sprintf('%.3f',round(dune.pairwise.adonis$F.Model,3))
  dune.pairwise.adonis$distance<-distance
  colnames(dune.pairwise.adonis)[1] <- "comparison"
  colnames(dune.pairwise.adonis)[4] <- "statistic (F)"
  colnames(dune.pairwise.adonis)[7] <- "pvalue"

  for (tt in 1:nrow(dune.pairwise.adonis)) {
    if (dune.pairwise.adonis$pvalue[tt] <= 0.001) dune.pairwise.adonis$sig[tt] <- '***'
    else if (dune.pairwise.adonis$pvalue[tt] <= 0.01) dune.pairwise.adonis$sig[tt] <- '**'
    else if (dune.pairwise.adonis$pvalue[tt] <= 0.05) dune.pairwise.adonis$sig[tt] <- '*'
    else dune.pairwise.adonis$sig[tt] <- ""
  }

  if (distance == 'bray') distance<-"Bray-curtis"
  if (distance == 'jaccard') distance<-"Jaccard"
  file2=paste(export_path,"/beta diversity/adonis/adonis_",distance,".txt",sep = "")
  write.table(dune.pairwise.adonis,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity based on Adonis_",distance," has been exported to","/",export_path,"",sep = "","\n")

  dune.pairwise.adonis<-subset(dune.pairwise.adonis, select=c(1,distance,4,7,8))
  dune.pairwise.adonis$pvalue<-sprintf('%.3f',round(dune.pairwise.adonis$pvalue,3))
  return(dune.pairwise.adonis)
}

"pairwise.mrpp" <-function(x,factors, sim.method, p.adjust.m){
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs=c()
  A = c()
  ObservedDelta=c()
  ExpectedDelta = c()
  p.value=c()
  for(elem in 1:ncol(co)){
    ad = mrpp(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),],
              factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))], permutations = 999,distance = sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    A = c(A,ad$A);
    ObservedDelta = c(ObservedDelta,ad$E.delta)
    ExpectedDelta = c(ExpectedDelta,ad$delta)
    p.value=c(p.value, ad$Pvalue)
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs, A, ObservedDelta, ExpectedDelta, p.value, p.adjusted)
  return(pairw.res)
}

"dune.pairwise.mrpp"<-function(otu,
                               p.adjust.m = "BH",
                               distance = 'bray',
                               export_path="cs2/microbial diversity analysis") {
  dir.create(paste(export_path,"/beta diversity/mrpp",sep = ""), recursive = TRUE)
  sim.method<-distance
  group<-group_generate(otu)
  otu<-t(otu)%>%data.frame()
  group$group<-factor(group$group,levels = unique(group$group))
  dune.pairwise.xmrpp<-pairwise.mrpp(otu, group$group, sim.method=sim.method, p.adjust.m= p.adjust.m)
  cat("\n","A: 'A > 0' states that the difference in between-group comparison is greater than that of within-group comparison.")
  cat("\n","Observed Delta: The smaller value of Observed Delta, the smaller within-group difference.")
  cat("\n","Expected Delta (∂): The smaller value of Expected Delta, the smaller between-group difference.")

  dune.pairwise.xmrpp$sig<-""
  for (tt in 1:nrow(dune.pairwise.xmrpp)) {
    if (dune.pairwise.xmrpp$p.adjusted[tt] <= 0.001) dune.pairwise.xmrpp$sig[tt] <- '***'
    else if (dune.pairwise.xmrpp$p.adjusted[tt] <= 0.01) dune.pairwise.xmrpp$sig[tt] <- '**'
    else if (dune.pairwise.xmrpp$p.adjusted[tt] <= 0.05) dune.pairwise.xmrpp$sig[tt] <- '*'
    else dune.pairwise.xmrpp$sig[tt] <- ""
  }

  dune.pairwise.xmrpp$ExpectedDelta<-sprintf('%.3f',round(dune.pairwise.xmrpp$ExpectedDelta,3))
  dune.pairwise.xmrpp$ObservedDelta<-sprintf('%.3f',round(dune.pairwise.xmrpp$ObservedDelta,3))
  dune.pairwise.xmrpp$A<-sprintf('%.3f',round(dune.pairwise.xmrpp$A,3))
  dune.pairwise.xmrpp$p.adjusted<-sprintf('%.3f',round(dune.pairwise.xmrpp$p.adjusted,3))


  dune.pairwise.xmrpp$distance<-distance
  colnames(dune.pairwise.xmrpp)[1] <- "comparison"
  colnames(dune.pairwise.xmrpp)[4] <- "statistic (∂)"
  colnames(dune.pairwise.xmrpp)[6] <- "pvalue"
  dune.pairwise.xmrpp<-subset(dune.pairwise.xmrpp,select=c(1,8,2,3,4,6,7))

  if (distance == 'bray') distance<-"Bray-curtis"
  if (distance == 'jaccard') distance<-"Jaccard"

  file2=paste(export_path,"/beta diversity/mrpp/mrpp_",distance,".txt",sep = "")
  write.table(dune.pairwise.xmrpp,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity based on MRPP_",distance," has been exported to","/",export_path,"",sep = "","\n")

  return(dune.pairwise.xmrpp)
}

"p.adjust.method" <- function() {
  ss<-c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
  st<-ss[1]
  for (tt in 2:length(ss)) {
    ts<-ss[tt]
    {
      st<-paste(st,ts)
    }
  }
  cat("\n","p.adjust.m (chosen):",st,sep=" ")
  message("\n"," p.adjust.m (chosen): ",st)
}

"pairwise.commdiff" <- function(submchat,
                                plot.data=FALSE,
                                p.adjust.m = "BH",
                                distance = 'bray',
                                export_path="cs2/microbial diversity analysis") {
  dir.create(paste(export_path,"/beta diversity",sep = ""), recursive = TRUE)
  suppressMessages(library(vegan))
  p.adjust.method()
  otu<-submchat$otu_table

  anosim.dat <- dune.pairwise.anosim(otu,
                                     p.adjust.m = p.adjust.m,
                                     distance = distance,
                                     export_path=export_path)
  adonis.dat <- dune.pairwise.adonis(otu,
                                     p.adjust.m = p.adjust.m,
                                     distance = distance,
                                     export_path=export_path)
  mrpp.dat <- dune.pairwise.mrpp(otu,
                                 p.adjust.m = p.adjust.m,
                                 distance = distance,
                                 export_path=export_path)

  comm.diff.dat<-cbind(anosim.dat[,c(1,3,4)],adonis.dat[,c(3,4)],mrpp.dat[,c(5,6)])
  colnames(comm.diff.dat)<-c("Comparison","r","P","F","P","∂","P")

  comm.diff.data<-rbind(colnames(comm.diff.dat),comm.diff.dat)
  colnames(comm.diff.data)[2:3]<-"ANOSIM"
  colnames(comm.diff.data)[4:5]<-"Adonis"
  colnames(comm.diff.data)[6:7]<-"MRPP"
  cat("\n","Three non-parametric dissimilarity analyses:")
  cat("\n","  Multiple response  permutation  procedure (MRPP)")
  cat("\n","  Analysis of similarity (ANOSIM)")
  cat("\n","  Permutational multivariate analysis of variance (Adonis&Permanova)")
  comm.diff.data[1,1]<-""
  #require(ggpubr)

  if (distance == 'bray') distance<-"Bray-curtis"
  if (distance == 'jaccard') distance<-"Jaccard"

  p<-ggpubr::ggtexttable(comm.diff.data,rows = NULL)

  subtitle<-paste0("Note: Three different permutation tests were performed (MRPP, ANOSIM and Adonis) on the basis of ",distance," distance.")%>%
    strwrap(width = 75) %>%
    paste(collapse = "\n")

  p<-p%>%
    ggpubr::tab_add_hline(at.row = c(1,3), row.side = "top",
                  linewidth = 3, linetype = 1) %>%
    ggpubr::tab_add_hline(at.row = c((length(comm.diff.data$Comparison)+1)),
                  row.side = "bottom", linewidth = 3, linetype = 1) %>%
    ggpubr::tab_add_vline(at.column = 2:ncol(comm.diff.data), column.side = "left",
                  from.row = 3, linetype = 2) %>%
    ggpubr::tab_add_title(text = "Significance tests of the networked communities between pairwise comparison",
                          family = "serif",
                          size = 10, face = "bold.italic")%>%
    ggpubr::tab_add_footnote(text = subtitle,
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic") %>%
    ggpubr::tab_add_footnote(text = "--Multiple response  permutation  procedure (MRPP)",
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic") %>%
    ggpubr::tab_add_footnote(text = "--Analysis of similarity (ANOSIM)",
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic") %>%
    ggpubr::tab_add_footnote(text = '--Permutational multivariate analysis of variance (Adonis or Permanova)',
                             family = "serif",padding = unit(1, "line"),
                             size = 10, face = "italic") %>%
    ggpubr::tab_add_footnote(text = "*All data processed by R were presented.",
                     family = "serif",padding = unit(1, "line"),
                     size = 10, face = "italic")

  print(p)



  file2=paste(export_path,"/beta diversity/Three permutation tests (",distance,").txt",sep = "")
  write.table(comm.diff.data,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","beta diversity based on three different permutation tests under ",distance," has been exported to","/",export_path,"",sep = "","\n")
  cat("\n","Copywriting (Title): Significance tests of the networked communities between pairwise comparison")
  cat("\n"," Copywriting (Note): Three different permutation tests were performed (MRPP, ANOSIM and Adonis) on the basis of ",distance," distance",sep = "")
  if (plot.data) comm.diff.data <- list(plot=p,data=comm.diff.data)
  return(comm.diff.data)
}

