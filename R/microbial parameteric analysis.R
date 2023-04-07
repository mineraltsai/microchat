"calcMicrochatParam" <- function(submchat,
                                 paramfile.select=c("gut_gene","gut_enzyme"),
                                 export_path="microbial parameteric analysis") {
  dir.create(export_path, recursive = TRUE)

  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  if (class(submchat$param_table)!="list") {
    param_table<-submchat$param_table
    paramfile.select<-NULL
    }
  if (class(submchat$param_table)=="list") {
    if (length(paramfile.select)==1) {
      param<-submchat$param_table
      param.select<-param[[which(names(param)==paramfile.select)]]
    } else {
      param<-submchat$param_table
      param.select<-param[[which(names(param)==paramfile.select[1])]]
      for (ps in paramfile.select[2:length(paramfile.select)]) {
        ps.int<-param[[which(names(param)==ps)]]
        {
          param.select<-cbind(param.select,ps.int)
        }
      }
    }
    param_table<-param.select
  }


  param_table<-tibble::rownames_to_column(param_table,var = "sample")
  param_table$group<-substr(param_table$sample,start = 1,stop = 2)
  param_table<-subset(param_table, select = c(which(colnames(param_table)=="group"),
                                              which(colnames(param_table)!="group")))


  param_table$group<-factor(param_table$group,levels = unique(param_table$group))

  ###normality test
  nor.data<-subset(param_table,select = -c(group,sample))
  nptest<-data.frame()
  for (i in 1:length(colnames(nor.data))) {
    sw<-norm.test(nor.data[,i])
    {
      nptest<-rbind(nptest,sw)
    }
  }

  cat("--------------------------------------------------------------------------------------")
  ###Homogeneity of variance
  variance.data<-subset(param_table,select = -c(sample))
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


  param_table$group<-factor(param_table$group,levels = unique(param_table$group))

  ###output statistics table
  all.index<-colnames(param_table)[3:length(colnames(param_table))]
  cat("\n Please use one of '",all.index,"' in the plot function !!!\n")

  summ_fcbv.all<-data.frame()
  for (index in all.index) {
    # calculate the mean and standard error of each group according to the alpha index
    summ_fcbv<- param_table %>%
      group_by(group) %>%
      plyr::summarise(
        mean = aggregate(as.formula(paste(index, " ~ group",sep = "")),data=param_table,FUN = "mean"),
        sd = aggregate(as.formula(paste(index, " ~ group",sep = "")),data=param_table,FUN = "sd")
      ) %>%data.frame()

    summ_fcbv1<-cbind(summ_fcbv$mean,summ_fcbv$sd)
    summ_fcbv1<-subset(summ_fcbv1, select=c(-3))
    colnames(summ_fcbv1)<-c("group","mean","sd")

    summ_fcbv1<-summ_fcbv1[1:length(unique(param_table$group)),]
    summ_fcbv1$group<-factor(summ_fcbv1$group,levels = unique(param_table$group))

    gpnum<-data.frame()
    for (gname in unique(param_table$group)) {
      gnum<-which(param_table$group==gname)%>%length()
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

  write.table(summ_fcbv.all,file = paste(export_path,"/parameters.stat.txt",sep = ""),row.names = FALSE,quote = FALSE, sep = "\t")
  message("\n","The param_table was used for statistical analysis. You could check it.","\n")

  message("\n","The param_table has been statistically analyzed. You could check it.","\n")

  microchatParamobj<-list(param_table,nptest,vttest,summ_fcbv.all,all.index,paramfile.select)

  names(microchatParamobj)[1]<-"param_table"
  names(microchatParamobj)[2]<-"np.test"
  names(microchatParamobj)[3]<-"vt.test"
  names(microchatParamobj)[4]<-"statistics"
  names(microchatParamobj)[5]<-"all.index"
  names(microchatParamobj)[6]<-"paramfile.select"

  class(microchatParamobj) <- c("microchat","list")

  return(microchatParamobj)
}





"calcMicrochatParamStat" <- function(microchatParamobj,
                                        select.index="ph",
                                        strictmod=FALSE,
                                        method=c("t.test","wilcox.test","anova","kruskal.test"),
                                        comparison=my_comparisons,
                                        export_path="microbial parameteric analysis") {


  method<-match.arg(method)
  dir.create(export_path, recursive = TRUE)

  if (class(microchatParamobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  alphadiv<-microchatParamobj$param_table
  alphadiv$group<-factor(alphadiv$group,levels = unique(alphadiv$group))
  ###output alpha table
  alphadiv.use<-reshape::melt(alphadiv)

  ### words match
  all.index<-microchatParamobj$all.index
  index<-all.index[match(tolower(select.index),tolower(all.index))]

  ###define new alpha table according to the selected alpha index
  data_poi<-alphadiv.use[which(alphadiv.use$variable==index),]
  data_poi$group <- factor(data_poi$group , levels = unique(data_poi$group))

  ###define new statistics table according to the selected alpha index
  data_err<-microchatParamobj$statistics
  data_err<-data_err[which(data_err$index==index),]
  data_err$group <- factor(data_err$group , levels = unique(data_err$group))


  ###assume the sample number varied
  groupnum<-unique(alphadiv$group)%>%length()
  samplenum<-unique(alphadiv$sample)%>%length()

  #my_comparisons<-calcComparisonOrder(my_comparisons,data_poi)

  if (!strictmod) {

    ##!strictmod
    ###assume the variance test of all data
    np.test<-microchatParamobj$np.test
    vp.test<-microchatParamobj$vt.test

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
    np.test<-microchatParamobj$np.test
    vp.test<-microchatParamobj$vt.test

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

  opc.data<-get.opc(microchatParamobj)
  #write.table(data.alpha,file = paste(export_path,"/alpha diversity.stat.diff.txt",sep = ""),row.names = FALSE,quote = FALSE, sep = "\t")
  message("The parameteric properities has been statistically analyzed. You could check it.")

  microchatParamStatobj<-list(data.alpha,data_poi,data_err,index,method,y.ratio,sig_label_new,opc.data)

  names(microchatParamStatobj)[1]<-"param.stats"
  names(microchatParamStatobj)[2]<-"data_poi"
  names(microchatParamStatobj)[3]<-"data_err"
  names(microchatParamStatobj)[4]<-"index"
  names(microchatParamStatobj)[5]<-"method"
  names(microchatParamStatobj)[6]<-"y.ratio"
  names(microchatParamStatobj)[7]<-"sig"
  names(microchatParamStatobj)[8]<-"opc.data"

  class(microchatParamStatobj) <- c("microchat","list")

  return(microchatParamStatobj)
}



"plotMicrochatParamBoxplot" <- function(microchatParamStatobj,
                                        xlabname=NULL,
                                        yaxis.italic=TRUE,
                                        errorbar.pos.adj=FALSE,
                                    errorbar.line.add=FALSE,
                                    seg=TRUE,
                                    add.spline=TRUE,
                                    errorbar.point.size=0.1,
                                    y.point.adj=NULL,
                                    color_group=colorCustom(5,pal = "gygn"),
                                    color_backgroud="grey90",
                                    export_path="microbial diversity analysis") {

  dir.create(export_path, recursive = TRUE)

  if (class(microchatParamStatobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data_poi<-microchatParamStatobj$data_poi
  data_err<-microchatParamStatobj$data_err
  data.alpha<-microchatParamStatobj$param.stats
  index<-microchatParamStatobj$index
  if (is.null(y.point.adj)) {
    y.ratio<-microchatParamStatobj$y.ratio
  } else {
    y.ratio<-microchatParamStatobj$y.ratio*(1+y.point.adj)
  }
  method<-microchatParamStatobj$method
  sig_label_new<-microchatParamStatobj$sig
  opc.data<-microchatParamStatobj$opc.data
  o.data<-opc.data[which(opc.data$index==index),]
  gnum<-ncol(opc.data.new)
  liner<-paste(colnames(o.data[gnum-3])," effect: ",o.data[gnum-3],sep = "")
  quad<-paste(colnames(o.data[gnum-2])," effect: ",o.data[gnum-2],sep = "")
  subtitile<-paste(liner,quad,sep = " ")

  colors<-color_group
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
                    aes(x = group, y = mean, group = index, ymin = mean-se, ymax = mean+se),
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

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


   if (yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="italic", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="italic", size=12))

    p<-p+
      theme(#panel.border = element_blank(),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_line(color='black'),
        panel.background = element_rect(fill = color_backgroud),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=10,family = "serif"),
        axis.text.x=element_text(colour = "black",size = 10,
                                 angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
        #strip.background =  element_blank(),
        legend.position = "none",aspect.ratio = 1)


  }

  if (method %in% c("anova","kruskal.test")){
    data.alpha$p.value<-data.alpha[,2]
    if (method=="anova") {
      gorder<-unique(data_poi$group)%>%as.character()
      sig_label_new$group<-ordered(sig_label_new$group,levels = gorder)
    p<-ggplot() +
      geom_errorbar(data = data_err,
                    aes(x = group, y = mean, group = index,
                        ymin = mean-se, ymax = mean+se),
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
      geom_text(data = sig_label_new,vjust=-0.5,
                aes(x = group,y = mean+se,label = alpha),
                size = 5,color = "black",family = "serif")

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


    if (yaxis.italic) p<-p+labs(y=index,
                                subtitle=subtitile,
                                title = paste(method,": p = ",
                                              keep.decmi(data.alpha$p.value),sep=""))+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(face = "bold.italic",family = "serif", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,
                                 subtitle=subtitile,
                                 title = paste(method,": p = ",
                                               keep.decmi(data.alpha$p.value),sep=""))+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
            title = element_text(face = "bold.italic",family = "serif", size=12))

    p<-p+
      theme(#panel.border = element_blank(),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_line(color='black'),
        panel.background = element_rect(fill = color_backgroud),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=10,family = "serif"),
        axis.text.x=element_text(colour = "black",size = 10,
                                 angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
        #strip.background =  element_blank(),
        legend.position = "none",aspect.ratio = 1)



    if (add.spline) p<-p+
      ggalt::geom_xspline(data= sig_label_new,size=0.5,color="grey50",
                                      aes(x=group, y=mean, group=index),spline_shape = 1)
    }

    if (method=="kruskal.test") {
      gorder<-unique(data_poi$group)%>%as.character()
      sig_label_new$group<-ordered(sig_label_new$group,levels = gorder)
      p<-ggplot(data=data_poi) +
        geom_errorbar(data = data_err,
                      aes(x = group, y = mean, group = index,
                          ymin = mean-se, ymax = mean+se),
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
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 5,color = "black",family = "serif")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")

      if (yaxis.italic) p<-p+labs(y=index,
                                  subtitle=subtitile,
                                  title = paste(method,": p = ",
                                                keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   subtitle=subtitile,
                                   title = paste(method,": p = ",
                                                 keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      p<-p+
        theme(#panel.border = element_blank(),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black'),
          panel.background = element_rect(fill = color_backgroud),
          #axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=10,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 10,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          #strip.background =  element_blank(),
          legend.position = "none",aspect.ratio = 1)


      if (add.spline) p<-p+
        ggalt::geom_xspline(data= sig_label_new,size=0.5,color="grey50",
                        aes(x=group, y=mean),spline_shape = 1)
      }

  }

  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-microchatParamStatobj$data_err$group%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p<-p+ scale_x_discrete(labels = xlabname)
  }

  p<-p+theme(title = element_text(size=8),aspect.ratio = 1)
  ggsave(paste(export_path,"/Parameter (",index,") .pdf",sep = ""),p)
  cat("Parametric properities barplot has been exported. Please check it.","\n")

  return(p)

}

"plotMicrochatParamMutiBoxplot" <- function(microchatParamobj,
                                            xlabname=NULL,
                                            yaxis.italic=TRUE,
                                            strictmod=TRUE,
                                            method="anova",
                                            comparison=my_comparisons,
                                            color_group=colorCustom(5,pal = "ywbu"),
                                            export_path="ss21/microbial parameteric analysis/liver_gene") {

  select.index<-microchatParamobj$all.index
  pp<-list()
  for (t in select.index) {
    microchatParamStatobj<-calcMicrochatParamStat(microchatParamobj,
                                                  select.index=t,
                                                  strictmod=strictmod,
                                                  method=method,
                                                  comparison=comparison,
                                                  export_path=export_path)

    p<-plotMicrochatParamBoxplot(microchatParamStatobj,
                                 xlabname=xlabname,
                                 yaxis.italic=yaxis.italic,
                                 errorbar.line.add=TRUE,
                                 errorbar.point.size=0.1,
                                 y.point.adj=0.1,
                                 seg=TRUE,
                                 add.spline=FALSE,
                                 color_group=color_group,
                                 color_backgroud="grey90",
                                 export_path=export_path)
    pp[[t]]<-p
  }

  library(patchwork)
  pm<-wrap_plots(pp)
  ggsave(paste(export_path,"/Muti-parameter",".pdf",sep = ""),pm)
  message("Muti-parametric properities barplot has been exported. Please check it.")

  return(pm)
}

"calcMicrochatParamTable" <- function(microchatParamobj,
                                      strictmod=TRUE,
                                      method="anova",
                                      comparison=my_comparisons,
                                      export_path="microbial parameteric analysis") {
  dir.create(export_path, recursive = TRUE)
  paramfile.select<-microchatParamobj$paramfile.select
  if (class(microchatParamobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  suppressMessages(library(reshape2))
  suppressMessages(library(microchat))

  data_err<-microchatParamobj$statistics
  data_err$summa<-paste(sprintf("%0.2f", data_err$mean),sprintf("%0.2f", data_err$se),sep = "±")

  param_tab<-subset(data_err,select=c(index,group,summa))%>%spread(key =group, value = summa )

  if (method %in% c("t.test","wilcox.test")) {
    select.index<-microchatParamobj$all.index
    sigtabx<-data.frame()
  for (t in select.index) {

    microchatParamStatobj<-calcMicrochatParamStat( microchatParamobj,
                                                  select.index=t,
                                                  strictmod=strictmod,
                                                  method=method,
                                                  comparison=comparison,
                                                  export_path=export_path)

    sigtab<-microchatParamStatobj$param.stats
    sigtab<-subset(sigtab,select=c(group1,group2,sig))
    sigtab$index<-t

    sigtab$group1<-factor(sigtab$group1,levels =unique(microchatParamStatobj$data_poi$group))
    sigtab$group2<-factor(sigtab$group2,levels =unique(microchatParamStatobj$data_poi$group))
    sigtab<-sigtab[order(sigtab$group1,sigtab$group2),]

    if (length(unique(sigtab$group1))==2) {
      gwrow<-which(as.character(unique(sigtab$group1)[2])==as.character(sigtab$group1))
      for (tk in gwrow) {
        if (sigtab$sig[tk]=="*") sigtab$sig[tk]<-"#"
        if (sigtab$sig[tk]=="**") sigtab$sig[tk]<-"#"
        if (sigtab$sig[tk]=="***") sigtab$sig[tk]<-"#"
        if (sigtab$sig[tk]=="****") sigtab$sig[tk]<-"#"
      }
    }

    sigtab$sig[which(sigtab$sig=="ns")]<-""

    {
      sigtabx<-rbind(sigtabx,sigtab)
    }
  }

    sigtabxx<-tidyr::spread(sigtabx,key=index,value=sig)
    sigtabxx<-subset(sigtabxx,select = -group1)
    sigtabxx<-t(sigtabxx)%>%data.frame()

    for (ak in 1:ncol(sigtabxx)) {
      group.select<-sigtabxx[1,][ak]%>%as.character()
      for (rrl in 1:nrow(param_tab)) {
        param_tab[rrl,which(group.select==colnames(param_tab))]<-paste(param_tab[rrl,which(group.select==colnames(param_tab))],
                                                                    sigtabxx[rrl+1,ak],sep = "")
      }
    }
  }

  ###variance analysis one-way
  if (method %in% c("anova","kruskal.test")) {
    select.index<-microchatParamobj$all.index
    sigtabx<-data.frame()
    for (t in select.index) {
      microchatParamStatobj<-calcMicrochatParamStat( microchatParamobj,
                                                     select.index=t,
                                                     strictmod=strictmod,
                                                     method=method,
                                                     comparison=comparison,
                                                     export_path=export_path)



      sigtab<-microchatParamStatobj$sig
      sigtab<-subset(sigtab,select=c(group,alpha))
      sigtab$index<-t
      sigtab$group<-factor(sigtab$group,levels =unique(microchatParamStatobj$data_poi$group))
      sigtab<-sigtab[order(sigtab$group),]
      {
        sigtabx<-rbind(sigtabx,sigtab)
      }
    }

    sigtabxx<-tidyr::spread(sigtabx,key=index,value=alpha)
    sigtabxx<-sigtabxx%>%tibble::column_to_rownames(var = "group")
    sigtabxx<-sigtabxx%>%t()%>%data.frame()
    param_tab<-param_tab%>%tibble::column_to_rownames(var = "index")

    rname<-rownames(param_tab)
    cname<-colnames(param_tab)
    for (rx in rname) {
      for (cx in cname) {
        param_tab[rx,cx]<-paste(param_tab[rx,cx],sigtabxx[rx,cx],sep = "")
      }
    }

    ###add opc prediction
    params_table<-microchatParamobj$param_table
    params_table<-params_table[,-1]
    params_table<-column_to_rownames(params_table,var = "sample")%>%data.frame()

    colname<-colnames(params_table)
    trt_id <-unique(colname)

    ttk<-sapply(trt_id,function(x){
      grep(x,colnames(params_table))
    })

    split_otu <-
      lapply(ttk,FUN = function(x){
        parax<-params_table[,x]%>%data.frame()
        rownames(parax)<-rownames(params_table)
        parax<-rownames_to_column(parax,var = "id")
        parax$group<-substr(parax$id,start = 1,stop = 2)
        parax$group<-factor(parax$group,levels = unique(parax$group))
        parax$id<-substr(parax$id,start = 3,stop=5)%>%as.numeric()
        colnames(parax)[2]<-"value"
        parax<-pivot_wider(parax,names_from = group,values_from = value)
        return(parax)
      })

    model.all<-c("Linear","Quadratic","Cubic","Quartic","Quintic","Hexic")
    fff<-lapply(lapply(split_otu, rma_opc),
                function(x){
                  opc_pre<-x$contrast_table
                  opc_pre<-subset(opc_pre,select=c(1,3,8))
                  colnames(opc_pre)<-c('Effect model',"Dose effect contribution","pvalue")
                  opc_pre$'Effect model'<-model.all[1:nrow(opc_pre)]
                  opc_pre$pvalue<-ifelse(opc_pre$pvalue<0.001,"<0.001",
                                         ifelse(opc_pre$pvalue<0.01,sprintf("%0.3f", round(opc_pre$pvalue,3)),
                                                ifelse(opc_pre$pvalue<0.05,sprintf("%0.3f", round(opc_pre$pvalue,3)),
                                                       sprintf("%0.3f", round(opc_pre$pvalue,3)))))
                  opc_pre$"Dose effect contribution"<-ifelse(opc_pre$"Dose effect contribution"<0.001,sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)),
                                                             ifelse(opc_pre$"Dose effect contribution"<0.01,sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)),
                                                                    ifelse(opc_pre$"Dose effect contribution"<0.05,sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)),
                                                                           sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)))))

                  opc_pre$comb<-paste(opc_pre$pvalue," (",opc_pre$"Dose effect contribution",") ",sep = "")
                  opc_pre<-subset(opc_pre,select=c(1,4))
                  return(opc_pre)
        })

    param_tab<-rownames_to_column(param_tab,var = "index")
    param_tab$index<-factor(param_tab$index,levels = names(fff))
    param_tab<-param_tab[order(param_tab$index),]
    param_tab[,(ncol(param_tab)+1):(ncol(param_tab)+nrow(fff[[1]]))]<-0
    addlen<-length(fff[[1]]$`Effect model`)
    colnames(param_tab)[(ncol(param_tab)-addlen+1):ncol(param_tab)]<-fff[[1]]$`Effect model`

    for (index in names(fff)) {
      param_tab[which(param_tab$index==index),(ncol(param_tab)-addlen+1):ncol(param_tab)]<-fff[[which(names(fff)==index)]]$comb
    }
  }



  if (!is.null(paramfile.select)) flie=paste(export_path,"/(",paramfile.select,")parameters_statistic_",method,".csv", sep = "")
  if (is.null(paramfile.select)) flie=paste(export_path,"/","parameters_statistic_",method,".csv", sep = "")
  write.csv(param_tab,file = flie,row.names = FALSE)
}




"plotMicrochatParamComplexBoxplot" <- function(submchat,
                                               strictmod=TRUE,
                                               method="t.test",
                                               comparison=my_comparisons,
                                               xlabname=NULL,
                                               yaxis.italic=TRUE,
                                               color_group=colorCustom(5,pal = "ywbu"),
                                               export_path="ss21/microbial parameteric analysis") {
  params<-submchat$param_table

  for (paramname in names(params)) {

    export_path1<-export_path
    export_path2<-paste(export_path1,"/",paramname,sep = "")

    microchatParamobj<-calcMicrochatParam(submchat,
                                          paramfile.select=paramname,
                                          export_path=export_path2)
    calcMicrochatParamTable(microchatParamobj,
                            strictmod=strictmod,
                            method=method,
                            comparison=comparison,
                            export_path=export_path1)

    plotMicrochatParamMutiBoxplot(microchatParamobj,
                                  xlabname=xlabname,
                                  yaxis.italic=yaxis.italic,
                                  strictmod=strictmod,
                                  method=method,
                                  comparison=comparison,
                                  color_group=color_group,
                                  export_path=export_path2)
  }
}



"plotMicrochatParamBarplot" <- function(microchatParamStatobj,
                                        xlabname=NULL,
                                        yaxis.italic=TRUE,
                                        errorbar.pos.adj=FALSE,
                                        errorbar.line.add=FALSE,
                                        seg=TRUE,
                                        add.spline=TRUE,
                                        errorbar.point.size=0.1,
                                        y.point.adj=NULL,
                                        color_group=colorCustom(5,pal = "gygn"),
                                        color_backgroud="grey90",
                                        export_path="microbial diversity analysis") {

  dir.create(export_path, recursive = TRUE)

  if (class(microchatParamStatobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data_poi<-microchatParamStatobj$data_poi
  data_err<-microchatParamStatobj$data_err
  data.alpha<-microchatParamStatobj$param.stats
  index<-microchatParamStatobj$index
  if (is.null(y.point.adj)) {
    y.ratio<-microchatParamStatobj$y.ratio
  } else {
    y.ratio<-microchatParamStatobj$y.ratio*(1+y.point.adj)
  }
  method<-microchatParamStatobj$method
  sig_label_new<-microchatParamStatobj$sig
  opc.data<-microchatParamStatobj$opc.data
  o.data<-opc.data[which(opc.data$index==index),]
  gnum<-ncol(opc.data.new)
  liner<-paste(colnames(o.data[gnum-3])," effect: ",o.data[gnum-3],sep = "")
  quad<-paste(colnames(o.data[gnum-2])," effect: ",o.data[gnum-2],sep = "")
  subtitile<-paste(liner,quad,sep = " ")

  colors<-color_group
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
                    aes(x = group, y = mean, group = index, ymin = mean-se, ymax = mean+se),
                    colour = "grey50",width=.25)+
      geom_bar(data = data_err,aes(x = group,y = mean,fill = group),
               stat = "identity",position = "dodge",
               width = 0.7,colour = "white")


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

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


    if (yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="bold.italic", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="bold.italic", size=12))

    p<-p+
      theme(#panel.border = element_blank(),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_line(color='black'),
        panel.background = element_rect(fill = color_backgroud),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=10,family = "serif"),
        axis.text.x=element_text(colour = "black",size = 10,
                                 angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
        #strip.background =  element_blank(),
        legend.position = "none",aspect.ratio = 1)


  }

  if (method %in% c("anova","kruskal.test")){
    data.alpha$p.value<-data.alpha[,2]
    if (method=="anova") {
      gorder<-unique(data_poi$group)%>%as.character()
      sig_label_new$group<-ordered(sig_label_new$group,levels = gorder)

      p<-ggplot() +
        geom_errorbar(data = data_err,
                      aes(x = group, y = mean, group = index,
                          ymin = mean-se, ymax = mean+se),
                      colour = "grey50",width=.25)+
        geom_bar(data = data_err,aes(x = group,y = mean,fill = group),
                 stat = "identity",position = "dodge",
                 width = 0.7,colour = "white") +
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 5,color = "black",family = "serif")


      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")


      if (yaxis.italic) p<-p+labs(y=index,
                                  subtitle=subtitile,
                                  title = paste(method,": p = ",
                                                keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   subtitle=subtitile,
                                   title = paste(method,": p = ",
                                                 keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      p<-p+
        theme(#panel.border = element_blank(),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black'),
          panel.background = element_rect(fill = color_backgroud),
          #axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=10,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 10,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          #strip.background =  element_blank(),
          legend.position = "none",aspect.ratio = 1)



      if (add.spline) p<-p+
        ggalt::geom_xspline(data= sig_label_new,size=0.5,color="grey50",
                            aes(x=group, y=mean, group=index),spline_shape = 1)
    }

    if (method=="kruskal.test") {
      gorder<-unique(data_poi$group)%>%as.character()
      sig_label_new$group<-ordered(sig_label_new$group,levels = gorder)

      p<-ggplot() +geom_errorbar(data = data_err,
                                 aes(x = group, y = mean, group = index,
                                     ymin = mean-se, ymax = mean+se),
                                 colour = "grey50",width=.25)+
        geom_bar(data = data_err,aes(x = group,y = mean,fill = group),
                 stat = "identity",position = "dodge",
                 width = 0.7,colour = "white") +

        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 5,color = "black",family = "serif")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")

      if (yaxis.italic) p<-p+labs(y=index,
                                  subtitle=subtitile,
                                  title = paste(method,": p = ",
                                                keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   subtitle=subtitile,
                                   title = paste(method,": p = ",
                                                 keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      p<-p+
        theme(#panel.border = element_blank(),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black'),
          panel.background = element_rect(fill = color_backgroud),
          #axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=10,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 10,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          #strip.background =  element_blank(),
          legend.position = "none",aspect.ratio = 1)


      if (add.spline) p<-p+
        ggalt::geom_xspline(data= sig_label_new,size=0.5,color="grey50",
                            aes(x=group, y=mean),spline_shape = 1)
    }

  }

  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-microchatParamStatobj$data_err$group%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p<-p+ scale_x_discrete(labels = xlabname)
  }
  p<-p+theme(title = element_text(size=8),aspect.ratio = 1)
  ggsave(paste(export_path,"/Parameter (",index,") .pdf",sep = ""),p)
  cat("Parametric properities barplot has been exported. Please check it.")

  return(p)

}


"plotMicrochatParamLineplot" <- function(microchatParamStatobj,
                                         xlabname=NULL,
                                         yaxis.italic=TRUE,
                                         errorbar.pos.adj=FALSE,
                                         errorbar.line.add=FALSE,
                                         seg=TRUE,
                                         add.spline=FALSE,
                                         errorbar.point.size=0.1,
                                         y.point.adj=NULL,
                                         color_group=colorCustom(5,pal = "gygn"),
                                         color_backgroud="grey90",
                                         export_path="microbial diversity analysis") {

  dir.create(export_path, recursive = TRUE)

  if (class(microchatParamStatobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data_poi<-microchatParamStatobj$data_poi
  data_err<-microchatParamStatobj$data_err
  data.alpha<-microchatParamStatobj$param.stats
  index<-microchatParamStatobj$index
  if (is.null(y.point.adj)) {
    y.ratio<-microchatParamStatobj$y.ratio
  } else {
    y.ratio<-microchatParamStatobj$y.ratio*(1+y.point.adj)
  }
  method<-microchatParamStatobj$method
  sig_label_new<-microchatParamStatobj$sig

  opc.data<-microchatParamStatobj$opc.data
  o.data<-opc.data[which(opc.data$index==index),]
  gnum<-ncol(opc.data.new)
  liner<-paste(colnames(o.data[gnum-3])," effect: ",o.data[gnum-3],sep = "")
  quad<-paste(colnames(o.data[gnum-2])," effect: ",o.data[gnum-2],sep = "")
  subtitile<-paste(liner,quad,sep = " ")
  colors<-color_group
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
              aes(x = group,y = value,fill = group))+
      geom_line(data= data_err,size=0.5,colour = "grey50",aes(x=group, y=mean, group=index)) +
      geom_errorbar(data = data_err,
                    aes(x = group, y = mean, group = index, ymin = mean-se, ymax = mean+se),
                    colour = "grey50",width=.25)


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

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


    if (yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="bold.italic", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="bold.italic", size=12))

    p<-p+
      theme(#panel.border = element_blank(),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_line(color='black'),
        panel.background = element_rect(fill = color_backgroud),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=10,family = "serif"),
        axis.text.x=element_text(colour = "black",size = 10,
                                 angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
        #strip.background =  element_blank(),
        legend.position = "none",aspect.ratio = 1)


  }

  if (method %in% c("anova","kruskal.test")){
    data.alpha$p.value<-data.alpha[,2]
    if (method=="anova") {
      gorder<-unique(data_poi$group)%>%as.character()
      sig_label_new$group<-factor(sig_label_new$group,levels = gorder)
      sig_label_new<-sig_label_new[order(sig_label_new$group),]

      p<-ggplot()+
        geom_line(data= data_err,size=0.5,colour = "grey50",aes(x=group, y=mean, group=index)) +
        geom_errorbar(data = data_err,
                      aes(x = group, y = mean, group = index,
                          ymin = mean-se, ymax = mean+se),
                      colour = "grey50",width=.25)+
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 5,color = "black",family = "serif")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")


      if (yaxis.italic) p<-p+labs(y=index,
                                  subtitle=subtitile,
                                  title = paste(method,": p = ",
                                                keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   subtitle=subtitile,
                                   title = paste(method,": p = ",
                                                 keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))


      p<-p+
        theme(#panel.border = element_blank(),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black'),
          panel.background = element_rect(fill = color_backgroud),
          #axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=10,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 10,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          #strip.background =  element_blank(),
          legend.position = "none",aspect.ratio = 1)



      if (add.spline) p<-p+
        ggalt::geom_xspline(data= sig_label_new,size=0.5,color="grey50",
                            aes(x=group, y=mean, group=index),spline_shape = 1)
    }

    if (method=="kruskal.test") {
      gorder<-unique(data_poi$group)%>%as.character()
      sig_label_new$group<-factor(sig_label_new$group,levels = gorder)
      sig_label_new<-sig_label_new[order(sig_label_new$group),]

      p<-ggplot() +geom_line(data= data_err,size=0.5,colour = "grey50",
                             aes(x=group, y=mean, group=index))+
        geom_errorbar(data = data_err,
                      aes(x = group, y = mean, group = index,
                          ymin = mean-se, ymax = mean+se),
                      colour = "grey50",width=.25)+
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 5,color = "black",family = "serif")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")
      if (yaxis.italic) p<-p+labs(y=index,
                                  subtitle=subtitile,
                                  title = paste(method,": p = ",
                                                keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   subtitle=subtitile,
                                   title = paste(method,": p = ",
                                                 keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))


      p<-p+
        theme(#panel.border = element_blank(),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black'),
          panel.background = element_rect(fill = color_backgroud),
          #axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          #axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=10,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 10,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          #strip.background =  element_blank(),
          legend.position = "none",aspect.ratio = 1)


      if (add.spline) p<-p+
        ggalt::geom_xspline(data= sig_label_new,size=0.5,color="grey50",
                            aes(x=group, y=mean),spline_shape = 1)
    }

  }

  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-microchatParamStatobj$data_err$group%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p<-p+ scale_x_discrete(labels = xlabname)
  }
  p<-p+theme(title = element_text(size=8),aspect.ratio = 1)
  ggsave(paste(export_path,"/Parameter (",index,") .pdf",sep = ""),p)
  cat("Parametric properities barplot has been exported. Please check it.")

  return(p)

}



"plotMicrochatParamMutiBarplot" <- function(microchatParamobj,
                                            xlabname=NULL,
                                            yaxis.italic=TRUE,
                                            strictmod=TRUE,
                                            method="anova",
                                            comparison=my_comparisons,
                                            color_group=colorCustom(5,pal = "ywbu"),
                                            export_path="ss21/microbial parameteric analysis/liver_gene") {

  select.index<-microchatParamobj$all.index
  pp<-list()
  for (t in select.index) {
    microchatParamStatobj<-calcMicrochatParamStat(microchatParamobj,
                                                  select.index=t,
                                                  strictmod=strictmod,
                                                  method=method,
                                                  comparison=comparison,
                                                  export_path=export_path)

    p<-plotMicrochatParamBarplot(microchatParamStatobj,
                                 xlabname=xlabname,
                                 yaxis.italic=yaxis.italic,
                                 errorbar.line.add=TRUE,
                                 errorbar.point.size=0.1,
                                 y.point.adj=0.1,
                                 seg=TRUE,
                                 add.spline=FALSE,
                                 color_group=color_group,
                                 color_backgroud="grey90",
                                 export_path=export_path)
    pp[[t]]<-p
  }

  library(patchwork)
  pm<-wrap_plots(pp)
  ggsave(paste(export_path,"/Muti-parameter",".pdf",sep = ""),pm)
  message("Muti-parametric properities barplot has been exported. Please check it.")

  return(pm)
}



"plotMicrochatParamComplexBarplot" <- function(submchat,
                                               strictmod=TRUE,
                                               method="t.test",
                                               comparison=my_comparisons,
                                               xlabname=NULL,
                                               yaxis.italic=TRUE,
                                               color_group=colorCustom(5,pal = "ywbu"),
                                               export_path="ss21/microbial parameteric analysis") {
  params<-submchat$param_table

  for (paramname in names(params)) {

    export_path1<-export_path
    export_path2<-paste(export_path1,"/",paramname,sep = "")

    microchatParamobj<-calcMicrochatParam(submchat,
                                          paramfile.select=paramname,
                                          export_path=export_path2)
    calcMicrochatParamTable(microchatParamobj,
                            strictmod=strictmod,
                            method=method,
                            comparison=comparison,
                            export_path=export_path1)

    plotMicrochatParamMutiBarplot(microchatParamobj,
                                  xlabname=xlabname,
                                  yaxis.italic=yaxis.italic,
                                  strictmod=strictmod,
                                  method=method,
                                  comparison=comparison,
                                  color_group=color_group,
                                  export_path=export_path2)
  }
}

"get.opc" <- function(microchatParamobj) {

  data_err<-microchatParamobj$statistics
  data_err$summa<-paste(sprintf("%0.2f", data_err$mean),sprintf("%0.2f", data_err$se),sep = "±")
  param_tab<-subset(data_err,select=c(index,group,summa))%>%spread(key =group, value = summa )

  ###add opc prediction
  params_table<-microchatParamobj$param_table
  params_table<-params_table[,-1]
  params_table<-column_to_rownames(params_table,var = "sample")%>%data.frame()

  colname<-colnames(params_table)
  trt_id <-unique(colname)
  ttk<-sapply(trt_id,function(x){
    grep(x,colnames(params_table))
  })

  split_otu <-
    lapply(ttk,FUN = function(x){
      parax<-params_table[,x]%>%data.frame()
      rownames(parax)<-rownames(params_table)
      parax<-rownames_to_column(parax,var = "id")
      parax$group<-substr(parax$id,start = 1,stop = 2)
      parax$group<-factor(parax$group,levels = unique(parax$group))
      parax$id<-substr(parax$id,start = 3,stop=5)%>%as.numeric()
      colnames(parax)[2]<-"value"
      parax<-pivot_wider(parax,names_from = group,values_from = value)
      return(parax)
    })

  model.all<-c("Linear","Quadratic","Cubic","Quartic","Quintic","Hexic")
  fff<-lapply(lapply(split_otu, rma_opc),
              function(x){
                opc_pre<-x$contrast_table
                opc_pre<-subset(opc_pre,select=c(1,3,8))
                colnames(opc_pre)<-c('Effect model',"Dose effect contribution","pvalue")
                opc_pre$'Effect model'<-model.all[1:nrow(opc_pre)]
                opc_pre$pvalue<-ifelse(opc_pre$pvalue<0.001,"<0.001",
                                       ifelse(opc_pre$pvalue<0.01,sprintf("%0.3f", round(opc_pre$pvalue,3)),
                                              ifelse(opc_pre$pvalue<0.05,sprintf("%0.3f", round(opc_pre$pvalue,3)),
                                                     sprintf("%0.3f", round(opc_pre$pvalue,3)))))
                opc_pre$"Dose effect contribution"<-ifelse(opc_pre$"Dose effect contribution"<0.001,sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)),
                                                           ifelse(opc_pre$"Dose effect contribution"<0.01,sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)),
                                                                  ifelse(opc_pre$"Dose effect contribution"<0.05,sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)),
                                                                         sprintf("%0.3f", round(opc_pre$"Dose effect contribution",4)))))

                opc_pre$comb<-paste(opc_pre$pvalue," (",opc_pre$"Dose effect contribution",") ",sep = "")
                opc_pre<-subset(opc_pre,select=c(1,4))
                return(opc_pre)
              })


  param_tab$index<-factor(param_tab$index,levels = names(fff))
  param_tab<-param_tab[order(param_tab$index),]
  param_tab[,(ncol(param_tab)+1):(ncol(param_tab)+nrow(fff[[1]]))]<-0
  addlen<-length(fff[[1]]$`Effect model`)
  colnames(param_tab)[(ncol(param_tab)-addlen+1):ncol(param_tab)]<-fff[[1]]$`Effect model`

  for (index in names(fff)) {
    param_tab[which(param_tab$index==index),(ncol(param_tab)-addlen+1):ncol(param_tab)]<-fff[[which(names(fff)==index)]]$comb
  }
  param_tab
  return(param_tab)
}
