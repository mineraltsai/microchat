"calcMicrochatParam" <- function(submchat,
                                 paramfile.select=c("gut_gene","gut_enzyme"),
                                 export_path="microbial parameteric analysis") {
  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  dir.create(export_path, recursive = TRUE)
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

  cat("--------------------------------------------------------------------------------------\n")
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
  #select.index<-all.index[2]
  index<-all.index[match(tolower(select.index),tolower(all.index))]

  ###define new alpha table according to the selected alpha index
  data_poi<-alphadiv.use[which(alphadiv.use$variable==index),]
  data_poi$group <- factor(data_poi$group , levels = unique(data_poi$group))

  data_max<-subset(data_poi,select=c(1,4))%>%
    group_by(group)%>%
    summarise_all(max)

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


      #if ( method== "t.test") {
        sig_label_new=NULL
        sta<-lapply(comparison, function(x){
          alphadiv.use.selected1<-data_poi[
            which(data_poi$group==x[1] | data_poi$group==x[2]),]

          norm.fit<-norm.test(alphadiv.use.selected1[,4])
          homn.fit<-variance.test(alphadiv.use.selected1,'value',alpha = 0.05)
          if (norm.fit$p.value>=0.05) {
            if (homn.fit$p.value>=0.05) {
              fit<-t.test(value~group, alphadiv.use.selected1, var.equal=TRUE,
                          paired = FALSE, alternative = 'two.sided')
            } else {
              fit<-t.test(value~group, alphadiv.use.selected1, var.equal=FALSE,
                          paired = FALSE, alternative = 'two.sided')
            }
            method="t.test"
          } else {
            fit<-wilcox.test(value~group, alphadiv.use.selected1,
                             exact = FALSE,correct=FALSE,
                             paired = FALSE, alternative = 'two.sided')
            method="wilcox.test"
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
     # }
      opc.data<-NULL
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
        fit <- aov(formula(paste(index, "~group")), data=variance.data)
        tukey_result <- agricolae::HSD.test(fit, trt="group", group=TRUE, console=FALSE)

        tukey <- tukey_result$groups
        #means <- tapply(input.data[,index], input.data$group, mean)
        #stds <- tapply(input.data[,index], input.data$group, sd)
        #means_str <- sprintf("%0.2f±%0.2f", means, stds)

        #
        tukey$group<-rownames(tukey)
        tukey$group<-factor(tukey$group,levels = levels(variance.data$group))
        tukey<-tukey[order(tukey$group),]
        tukey[,index]<-tukey$groups
        tukey<-subset(tukey,select = -groups)

        tukey_new<-merge(tukey,data_max,by="group")
        colnames(tukey_new)[2]<-"alpha"
        tukey_new$group <- factor(tukey_new$group , levels = unique(tukey$group))
        tukey_new<-merge(tukey_new,data_err,by="group")
        tukey_new$valuey<-tukey_new$value+tukey_new$se
        sig_label_new<-tukey_new

        aov_table <- anova(fit)
        pvalue<-  aov_table$`Pr(>F)`[1]
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


        comparison<-with(variance.data,agricolae::kruskal(eval(parse(text=index)),group,group=TRUE))

        colnames(comparison$means)[1]<-"mean"
        colnames(comparison$groups)[1]<-index

        input.data<-comparison$means
        input.data<-input.data%>%tibble::rownames_to_column("group")
        input.data<-subset(input.data, select = c(1,4,2))

        input.let<-comparison$groups
        input.let<-input.let%>%tibble::rownames_to_column("group")
        input.data<-merge(input.data,input.let,by="group")
        input.data<-input.data[,-4]
        colnames(input.data)[4]<-"Letters"
        input.data[,"variable"]<-index

        input.data<-merge(input.data,data_max,by="group")
        sig_label_new<-input.data
        samplenum<-variance.data$group%>%table()%>%max()
        sig_label_new$valuey<-sig_label_new$mean+sig_label_new$std*sqrt(samplenum)
        sig_label_new$alpha<-sig_label_new$Letters


        aovtest=comparison$statistics$p.chisq
        if (aovtest<0.001) aovtest=0.001
        data.alpha<-data.frame(
          index=index,
          p.value=aovtest,
          method=method
        )
        y.ratio=NULL
      }
      opc.data<-get.opc(microchatParamobj)
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

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

      #if (index %in% match.index) {
      #  method= "t.test"
      #} else {
      #  method= "wilcox.test"
      #}

      #if ( method== "t.test") {
        sig_label_new=NULL
        sta<-lapply(comparison, function(x){
          alphadiv.use.selected1<-data_poi[
            which(data_poi$group==x[1] | data_poi$group==x[2]),]

          norm.fit<-norm.test(alphadiv.use.selected1[,4])
          if (norm.fit$p.value>=0.05) {
            homn.fit<-variance.test(alphadiv.use.selected1,'value',alpha = 0.05)
            if (homn.fit$p.value>=0.05) {
              fit<-t.test(value~group, alphadiv.use.selected1, var.equal=TRUE,
                          paired = FALSE, alternative = 'two.sided')
            } else {
              fit<-t.test(value~group, alphadiv.use.selected1, var.equal=FALSE,
                          paired = FALSE, alternative = 'two.sided')
            }
            method="t.test"
          } else {
            #homn.fit<-fligner.test(value ~ group, data = alphadiv.use.selected1)
            fit<-wilcox.test(value~group, alphadiv.use.selected1,
                             exact = FALSE,correct=FALSE,
                             paired = FALSE, alternative = 'two.sided')
            method="wilcox.test"
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
      #}

      opc.data<-NULL
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
        fit <- aov(formula(paste(index, "~group")), data=variance.data)
        tukey_result <- agricolae::HSD.test(fit, trt="group", group=TRUE, console=FALSE)

        tukey <- tukey_result$groups
        #means <- tapply(input.data[,index], input.data$group, mean)
        #stds <- tapply(input.data[,index], input.data$group, sd)
        #means_str <- sprintf("%0.2f±%0.2f", means, stds)

        #
        tukey$group<-rownames(tukey)
        tukey$group<-factor(tukey$group,levels = levels(variance.data$group))
        tukey<-tukey[order(tukey$group),]
        tukey[,index]<-tukey$groups
        tukey<-subset(tukey,select = -groups)

        tukey_new<-merge(tukey,data_max,by="group")
        colnames(tukey_new)[2]<-"alpha"
        tukey_new$group <- factor(tukey_new$group , levels = unique(tukey$group))
        tukey_new<-merge(tukey_new,data_err,by="group")
        tukey_new$valuey<-tukey_new$value+tukey_new$se
        sig_label_new<-tukey_new

        aov_table <- anova(fit)
        pvalue<-  aov_table$`Pr(>F)`[1]
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


        comparison<-with(variance.data,agricolae::kruskal(eval(parse(text=index)),group,group=TRUE))

        colnames(comparison$means)[1]<-"mean"
        colnames(comparison$groups)[1]<-index

        input.data<-comparison$means
        input.data<-input.data%>%tibble::rownames_to_column("group")
        input.data<-subset(input.data, select = c(1,4,2))

        input.let<-comparison$groups
        input.let<-input.let%>%tibble::rownames_to_column("group")
        input.data<-merge(input.data,input.let,by="group")
        input.data<-input.data[,-4]
        colnames(input.data)[4]<-"Letters"
        input.data[,"variable"]<-index

        input.data<-merge(input.data,data_max,by="group")
        sig_label_new<-input.data
        samplenum<-variance.data$group%>%table()%>%max()
        sig_label_new$valuey<-sig_label_new$mean+sig_label_new$std*sqrt(samplenum)
        sig_label_new$alpha<-sig_label_new$Letters


        aovtest=comparison$statistics$p.chisq
        if (aovtest<0.001) aovtest=0.001
        data.alpha<-data.frame(
          index=index,
          p.value=aovtest,
          method=method
        )
        y.ratio=NULL
      }
      opc.data<-get.opc(microchatParamobj)
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

    }

  }


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

"linear_equa_fitTest" <-function(data_err,start_pos=1){
  asd<-data_err
  asd$gg<-start_pos:(nrow(data_err)+start_pos-1)
  asd$gg<-asd$gg%>%as.numeric()
  # 拟合二次函数
  model <- lm(asd$mean ~ poly(asd$gg, 1, raw = TRUE), data = asd)

  ###5.非线性回归的方程
  mylr = function(x,y){
    ###绘点图
    #plot(x,y)

    x_mean = mean(x)
    y_mean = mean(y)

    a = coef(model)[2]
    b = coef(model)[1]

    # 构造非线性回归方程
    f <- a *x +b

    ###点图按照上述方程拟合曲线
    #curve(c + b*x + a*x^2, add = TRUE)
    ###生成总平方和
    sst = sum((y-y_mean)^2)
    ##残差平方和
    sse = sum((y-f)^2)
    ##回归平方和
    ssr = sum((f-y_mean)^2)

    result = c(a,b,sst,sse,ssr)
    names(result) = c('a','b','sst','sse','ssr')
    return(result)
  }

  ##横坐标为时间
  x = asd$gg
  ##纵坐标为表达量
  y = asd$mean
  ###代入公式
  f = mylr(x,y)

  ###计算拟合度，回归平方和除以总平方和
  R2= f['ssr']/f['sst']

  f = mylr(x,y)
  a<-f['a']
  b<-f['b']

  asd$y_new<-a*asd$gg+b

  kk<-list(R2=R2,a=a,b=b,formula=f,asd=asd)

  return(kk)
}
"obtainSlope" <-function(data_err){
  asd<-data_err
  asd$gg<-1:nrow(data_err)
  asd$gg<-asd$gg%>%as.numeric()
  # 拟合二次函数
  model <- lm(asd$mean ~ poly(asd$gg, 1, raw = TRUE), data = asd)

  ###5.非线性回归的方程
  mylr = function(x,y){
    ###绘点图
    #plot(x,y)

    x_mean = mean(x)
    y_mean = mean(y)

    a = coef(model)[2]
    b = coef(model)[1]

    # 构造非线性回归方程
    f <- a *x +b

    ###点图按照上述方程拟合曲线
    #curve(c + b*x + a*x^2, add = TRUE)
    ###生成总平方和
    sst = sum((y-y_mean)^2)
    ##残差平方和
    sse = sum((y-f)^2)
    ##回归平方和
    ssr = sum((f-y_mean)^2)

    result = c(a,b,sst,sse,ssr)
    names(result) = c('a','b','sst','sse','ssr')
    return(result)
  }

  ##横坐标为时间
  x = asd$gg
  ##纵坐标为表达量
  y = asd$mean
  ###代入公式
  f = mylr(x,y)
  a<-f['a']

  return(a)
}
"getCrossCoord" <- function(data1,data2) {
  data_err_test1<-data1
  data_err_test2<-data2
  a1 <- data_err_test1$a
  b1 <- data_err_test1$b
  a2 <- data_err_test2$a
  b2 <- data_err_test2$b

  x <- (b2 - b1) / (a1 - a2)
  y <- a1 * x + b1
  return(list(x=x,y=y))
}
"qua_equa_fitTest" <-function(data_err){
  asd<-data_err
  asd$gg<-1:nrow(data_err)
  asd$gg<-asd$gg%>%as.numeric()
  # 拟合二次函数
  model <- lm(asd$mean ~ poly(asd$gg, 2, raw = TRUE), data = asd)

  ###5.非线性回归的方程
  mylr = function(x,y){
    ###绘点图
    #plot(x,y)

    x_mean = mean(x)
    y_mean = mean(y)

    a = coef(model)[3]
    b = coef(model)[2]
    c = coef(model)[1]

    # 构造非线性回归方程
    f <- a * x^2 + b *x +c

    ###点图按照上述方程拟合曲线
    #curve(c + b*x + a*x^2, add = TRUE)
    ###生成总平方和
    sst = sum((y-y_mean)^2)
    ##残差平方和
    sse = sum((y-f)^2)
    ##回归平方和
    ssr = sum((f-y_mean)^2)

    result = c(a,b,c,sst,sse,ssr)
    names(result) = c('a','b','c','sst','sse','ssr')
    return(result)
  }

  ##横坐标为时间
  x = asd$gg
  ##纵坐标为表达量
  y = asd$mean
  ###代入公式
  f = mylr(x,y)

  ###计算拟合度，回归平方和除以总平方和
  R2= f['ssr']/f['sst']

  opt.dose<-(-coef(model)[2])/(2*coef(model)[3])
  opt.dose

  f = mylr(x,y)
  a<-f['a']
  b<-f['b']
  c<-f['c']

  asd$y_new<-a*asd$gg^2+b*asd$gg+c

  kk<-list(R2=R2,a=a,b=b,c=c,
           opt.dose=opt.dose,formula=f,asd=asd)

  return(kk)
}

"cub_equa_fitTest" <-function(data_err){
  asd<-data_err
  asd$gg<-1:nrow(data_err)
  asd$gg<-asd$gg%>%as.numeric()
  # 拟合二次函数
  model <- lm(asd$mean ~ poly(asd$gg, 3, raw = TRUE), data = asd)

  # 提取多项式系数
  coefficients <- coef(model)

  # 定义三次函数
  func <- function(x) {
    return(coefficients[1] + coefficients[2]*x + coefficients[3]*x^2 + coefficients[4]*x^3)
  }

  ###5.非线性回归的方程
  mylr = function(x,y){
    ###绘点图
    #plot(x,y)
    #line(x,y)
    x_mean = mean(x)
    y_mean = mean(y)

    a = coef(model)[3]
    b = coef(model)[2]
    c = coef(model)[1]
    d = coef(model)[4]

    # 构造非线性回归方程
    f <- d * x^3 + a * x^2 + b *x +c

    ###点图按照上述方程拟合曲线
    #curve(c + b*x + a*x^2, add = TRUE)
    ###生成总平方和
    sst = sum((y-y_mean)^2)
    ##残差平方和
    sse = sum((y-f)^2)
    ##回归平方和
    ssr = sum((f-y_mean)^2)

    result = c(d,a,b,c,sst,sse,ssr)
    names(result) = c('d','a','b','c','sst','sse','ssr')
    return(result)
  }

  ##横坐标为时间
  x = asd$gg
  ##纵坐标为表达量
  y = asd$mean
  ###代入公式
  f = mylr(x,y)

  ###计算拟合度，回归平方和除以总平方和
  R2= f['ssr']/f['sst']

  opt.dose<-optimize(func, interval = c(1, nrow(data_err)),maximum = TRUE)$maximum

  f = mylr(x,y)
  d<-f['d']
  a<-f['a']
  b<-f['b']
  c<-f['c']

  asd$y_new<- d*asd$gg^3+a*asd$gg^2+b*asd$gg+c


  # 预测拟合的值
  x_pred <- seq(1,nrow(data_err), length.out = 100)
  y_pred <- d*x_pred^3+a*x_pred^2+b*x_pred+c

  kk<-list(R2=R2,d=d,a=a,b=b,c=c,x_pred=x_pred,y_pred=y_pred,
           opt.dose=opt.dose,formula=f,asd=asd)

  return(kk)
}

"plotMicrochatParamBoxplot" <- function(microchatParamStatobj,
                                        geom.line=TRUE,
                                        geom.line.size=0.5,
                                        ylim.fold=1.1,
                                        xlabname=NULL,
                                        yaxis.italic=TRUE,
                                        errorbar.pos.adj=TRUE,
                                    errorbar.line.add=FALSE,
                                    seg=TRUE,
                                    panel.border=c("none","normal","all","axis"),
                                    errorbar.point.size=0.1,
                                    y.point.adj=NULL,
                                    color_group=colorCustom(5,pal = "gygn"),
                                    color_backgroud=NA,
                                    axis.x.angle.adjust=FALSE,
                                    mytheme=NULL,
                                    spline.params=NULL,
                                    export_path="microbial diversity analysis") {
  panel.border<-match.arg(panel.border)
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

  if (method %in% c("anova","kruskal.test")) {
  opc.data <- microchatParamStatobj$opc.data
  o.data <- opc.data[which(opc.data$index == index), ]
  gnum <- ncol(o.data)
  gnum/2
  liner <- paste(colnames(o.data[gnum/2 + 2]), " effect: ", # " effect: "
                 o.data[gnum/2 + 2], sep = "")
  quad <- paste(colnames(o.data[gnum/2 +3]), " effect: ",
                o.data[gnum/2 +  3], sep = "")
  subtitile <- paste(liner, "\n",quad, sep = "")
  }

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

    p<-ggplot(data=data_poi,aes(x=group, y=value))
    if (geom.line) p<-p+geom_line(data= data_err,size=geom.line.size,colour = "black",
                                  aes(x=group, y=mean, group=index))
    p<-p+geom_errorbar(data = data_err,
                       aes(x = group, y = mean, group = index,
                           ymin = mean-sd, ymax = mean+sd),
                       colour = "grey50",width=.25)+
      geom_boxplot(aes(fill=group),
        outlier.shape = NA,
                   width = 0.5,
                   color = "white")+
      geom_jitter(data=data_poi,size=3,alpha=0.5,
                  aes(group=index,color=group))+
      geom_point(aes(group=index,color=group))


    if (!seg) data.seg <- calcSegmentOrder(data.alpha, y.ratio,
                                           data_poi)
    if (seg) data.seg <-calcorderx(data.alpha, y.ratio,
                                   data_poi)

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

      data.seg$y1[tt]<-y1.use
      data.seg$y2[tt]<-y2.use
      data.seg$y3[tt]<-y3.use
      data.seg$y4[tt]<-y4.use

    }

    for (tt in 1:nrow(data.seg)) {
      x1<-data.seg$order1[tt]
      x2<-data.seg$order2[tt]

      y1.use<-data.seg$y1[tt]
      y2.use<-data.seg$y2[tt]
      y3.use<-data.seg$y3[tt]
      y4.use<-data.seg$y4[tt]

      label.use=data.seg$sig[tt]

      p<-p+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x1, y = y1.use, ymin = y1.use, ymax = y1.use)+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x2, y = y2.use, ymin = y2.use, ymax = y2.use)+
        annotate("segment",color="grey50",
                 x = x1, y = y1.use,
                 xend = x1, yend = y3.use)+
        annotate("segment",color="grey50",
                 x = x2, y = y2.use,
                 xend = x2, yend = y3.use)+
        annotate("segment",color="grey50",
                 x = x1, y = y3.use,
                 xend = x2, yend = y3.use)+
        annotate("text",family="serif",x = (x1+x2)/2, y = y3.use+y.ratio*2/3,label=label.use)

      if (errorbar.line.add) p<-p+annotate("segment", color="grey50",x = x1-0.05, y = y1.use, xend = x1+0.05, yend = y1.use)+
        annotate("segment",color="grey50", x = x2-0.05, y = y2.use, xend = x2+0.05, yend = y2.use)

    }

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


   if (yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="italic", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

      p<-ggplot(data=data_poi,aes(x=group, y=value))
      if (geom.line) p<-p+geom_line(data= data_err,size=geom.line.size,colour = "black",
                                    aes(x=group, y=mean, group=index))
      p<-p+geom_errorbar(data = data_err,
                         aes(x = group, y = mean, group = index,
                             ymin = mean-sd, ymax = mean+sd),
                         colour = "grey50",width=.25)+
      geom_boxplot(data=data_poi,
                   aes(x = group,y = value,fill = group),
                   outlier.shape = NA,
                   width = 0.5,
                   color = "white")+
      geom_jitter(data=data_poi,size=2.5,alpha=0.75,
                  aes(x=group, y=value, group=index,color=group))+
      geom_point(data=data_poi,
                 aes(x=group, y=value, group=index,color=group)) +
      geom_text(data = sig_label_new,#vjust=-0.5,
                nudge_y=1/20*max(data_poi$value),
                aes(x = group,
                    y = value,  ## mean+se
                    label = alpha),
                size = 3,color = "black",family = "serif", fontface="bold")

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


    if (yaxis.italic) p<-p+labs(y=index,
                                title = paste(
                                  paste(method,": p = ",
                                        keep.decmi(data.alpha$p.value),sep=""),"\n",
                                  subtitile,sep = ""
                                ))+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(face = "bold.italic",family = "serif", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,
                                 title = paste(
                                   paste(method,": p = ",
                                         keep.decmi(data.alpha$p.value),sep=""),"\n",
                                   subtitile,sep = ""
                                 ))+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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

    }

    if (method=="kruskal.test") {
      gorder<-unique(data_poi$group)%>%as.character()
      sig_label_new$group<-ordered(sig_label_new$group,levels = gorder)
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""


      p<-ggplot(data=data_poi,aes(x=group, y=value))
      if (geom.line) p<-p+geom_line(data= data_err,size=geom.line.size,colour = "black",
                  aes(x=group, y=mean, group=index))
      p<-p+geom_errorbar(data = data_err,
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
        geom_text(data = sig_label_new,#vjust=-0.5,
                  nudge_y=1/20*max(data_poi$value),
                  aes(x = group,
                      y = value, #### mean+se
                      label = alpha),
                  size = 3,color = "black",family = "serif", fontface="bold")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")

      if (yaxis.italic) p<-p+labs(y=index,
                                  title = paste(
                                    paste(method,": p = ",
                                          keep.decmi(data.alpha$p.value),sep=""),"\n",
                                    subtitile,sep = ""
                                  ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   title = paste(
                                     paste(method,": p = ",
                                           keep.decmi(data.alpha$p.value),sep=""),"\n",
                                     subtitile,sep = ""
                                   ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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

      }
  }

  ##add fitted curve for  quadratic model
  if (method %in% c("anova","kruskal.test") & !is.null(spline.params)){
    opc.thres=spline.params$opc.thres
    spline.size = spline.params$spline.size
    spline.color=spline.params$spline.color
    spline.shape =spline.params$spline.shape
    spline.type=spline.params$spline.type
    opt.interval=spline.params$opt.interval
    optimal.pt.size=spline.params$optimal.pt.size
    optimal.pt.color=spline.params$optimal.pt.color
    optimal.pt.type=spline.params$optimal.pt.type

    vline.color=spline.params$vline.color
    vline.size=spline.params$vline.size
    vline.type=spline.params$vline.type

    csal<-strsplit(o.data$Linear,"( )")[[1]][2]
    lin.num<-substr(csal,start = 2,stop =nchar(csal)-1 )%>%as.numeric()

    csa<-strsplit(o.data$Quadratic,"( )")[[1]][2]
    qua.num<-substr(csa,start = 2,stop =nchar(csa)-1 )%>%as.numeric()

    if (!is.na(lin.num)) {

    if ("Cubic" %in% colnames(o.data)) {
      csac<-strsplit(o.data$Cubic,"( )")[[1]][2]
      cub.num<-substr(csac,start = 2,stop =nchar(csac)-1 )%>%as.numeric()

      all.opc.r<-c(lin.num,qua.num,cub.num)

    } else {
      all.opc.r<-c(lin.num,qua.num)
    }

    #拟合三次函数
      if ("Cubic" %in% colnames(o.data)) if (cub.num>opc.thres & max(all.opc.r)==cub.num) {
      data_err_test<-cub_equa_fitTest(microchatParamStatobj$data_err)

      asd<-data_err_test$asd
      opt.dose<-data_err_test$opt.dose
      a<-data_err_test$a
      b<-data_err_test$b
      c<-data_err_test$c
      d<-data_err_test$d
      opt.inclusion<-opt.interval[as.integer(opt.dose)]+(opt.dose-as.integer(opt.dose))*(opt.interval[as.integer(opt.dose)+1]-opt.interval[as.integer(opt.dose)])

      l <- list(d = as.numeric(format(d, digits = 2)),
                a = as.numeric(format(a, digits = 2)),
                b = as.numeric(format(b, digits = 2)),
                c = as.numeric(format(c, digits = 2)),
                dose = as.numeric(format(opt.inclusion, digits = 2)),
                r2 = as.numeric(format(data_err_test$R2, digits = 2)))
      #eq<-substitute(italic(y) ==  d %.% (italic(x)-dose)^3 + a %.% (italic(x)-dose)^2 + b %.% (italic(x)-dose) + c~","~italic(R)^2~"="~r2,l)
      eq<-""

      p<-p + geom_line(data = data.frame(x = data_err_test$x_pred,
                                         y = data_err_test$y_pred),
                       aes(x, y),
                       color = spline.color,
                       linewidth = spline.size/2,
                       linetype=spline.type,
                       lineend = "round")+
        annotate(geom="segment",x = opt.dose, xend = opt.dose, size=vline.size,linetype=vline.type,
                 lineend = "round",
                 y =-Inf, yend = d*opt.dose^3+ a*opt.dose^2+b*opt.dose+c,color=vline.color)+
        geom_point(aes(x=opt.dose,y=d*opt.dose^3+a*opt.dose^2+b*opt.dose+c),
                   shape=optimal.pt.type,
                   size=optimal.pt.size,color=optimal.pt.color)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(d*opt.dose^3+a*opt.dose^2+b*opt.dose+c)),
                                 aes(x=opt.dose,y=y,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                 arrow = arrow(length = unit(0.01, "npc")),family="serif",size=2.5,
                                 box.padding = 1)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)*0.8),
                                 family="serif",size=2,
                                 aes(x = opt.dose, y = -Inf, label = as.character(as.expression(eq))),parse = TRUE)
    }
    #拟合二次函数
    if (qua.num>opc.thres & max(all.opc.r)==qua.num) {
      data_err_test<-qua_equa_fitTest(microchatParamStatobj$data_err)
      asd<-data_err_test$asd
      opt.dose<-data_err_test$opt.dose
      a<-data_err_test$a
      b<-data_err_test$b
      c<-data_err_test$c
      opt.inclusion<-opt.interval[as.integer(opt.dose)]+(opt.dose-as.integer(opt.dose))*(opt.interval[as.integer(opt.dose)+1]-opt.interval[as.integer(opt.dose)])

      l <- list(a = as.numeric(format(a, digits = 3)),
                b = as.numeric(format(b, digits = 3)),
                c=as.numeric(format(c, digits = 3)),
                dose = as.numeric(format(opt.inclusion, digits = 2)),
                r2 = as.numeric(format(data_err_test$R2, digits = 3)))
      #eq<-substitute(italic(y) == a %.% (italic(x)-dose)^2 + b %.% (italic(x)-dose) + c~","~italic(R)^2~"="~r2,l)
      eq<-""

      p<-p + ggalt::geom_xspline(data = asd,
                                 linetype=spline.type,lineend = "round",
                                 size = spline.size, color = spline.color,
                                 aes(x = gg, y = y_new), spline_shape = spline.shape) +
        annotate(geom="segment",x = opt.dose, xend = opt.dose, size=vline.size,linetype=vline.type,
                 lineend = "round",
                 y =-Inf, yend =  a*opt.dose^2+b*opt.dose+c,color=vline.color)+
        geom_point(aes(x=opt.dose,y=a*opt.dose^2+b*opt.dose+c),
                   shape=optimal.pt.type,
                   size=optimal.pt.size,color=optimal.pt.color)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)),
                                 aes(x=opt.dose,y=y,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                 arrow = arrow(length = unit(0.01, "npc")),family="serif",size=2.5,
                                 box.padding = 1)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)*0.8),
                                 family="serif",size=2,
                                 aes(x = opt.dose, y = -Inf, label = as.character(as.expression(eq))),parse = TRUE)
    }
    #拟合线性模型
    if (lin.num>opc.thres & max(all.opc.r)==lin.num) {
      data_err_test<-linear_equa_fitTest(microchatParamStatobj$data_err)

      if (data_err_test$a>=0) {
        message("The fitted curve was drawn as user required. However, the slope of the fitted curve is higher than 0, so the fitted curve will not display.")
      } else {
        message("The fitted curve was drawn as user required. The slope of the fitted curve is lower than 0, so the fitted curve will be divided into two segments.")
        message("Namely, fitted curve for broken-line model is showing  !!!")

        slope<-c()
        for (t in 1:(nrow(data_err_test$asd)-1)) {
          slope<-c(slope,obtainSlope(data<-data_err_test$asd[t:(t+1),]))
        }
        slope_a<-slope[which(slope<0)]
        slope_a<-slope_a[which(slope_a==min(slope_a))]
        spot_start<-which(slope==slope_a)
        if (which(data_err_test$asd$mean==max(data_err_test$asd$mean))!=as.numeric(spot_start)) {
          spot_start<-which(data_err_test$asd$mean==max(data_err_test$asd$mean))
        } else {
          spot_start<-spot_start
        }
        if (spot_start ==2 | spot_start==3) {
          spot_end<-spot_start+1

          data_err_test1<-linear_equa_fitTest(microchatParamStatobj$data_err[1:spot_start,])
          data_err_test2<-linear_equa_fitTest(microchatParamStatobj$data_err[spot_end:nrow(microchatParamStatobj$data_err),],spot_end)

          crossdot<-getCrossCoord(data_err_test1,data_err_test2)

          if (crossdot$x>1 & crossdot$x <nrow(data_err_test$asd)) {
            asd1<-data_err_test1$asd
            a1<-data_err_test1$a
            b1<-data_err_test1$b
            asd2<-data_err_test2$asd
            asd2$gg<-(max(asd1$gg)+1):nrow(data_err_test$asd)
            a2<-data_err_test2$a
            b2<-data_err_test2$b

            asd1_a<-asd1[1:2,]
            asd1_a$gg[1]<-min(asd1$gg)-0.3
            asd1_a$y_new[1]<-a1*asd1_a$gg[1]+b1
            asd1_a$gg[2]<-crossdot$x+0.3
            asd1_a$y_new[2]<-a1*asd1_a$gg[2]+b1
            asd1<-rbind(asd1,asd1_a)

            asd2_a<-asd2[1:2,]
            asd2_a$gg[1]<-max(asd2$gg)+0.3
            asd2_a$y_new[1]<-a2*asd2_a$gg[1]+b2
            asd2_a$gg[2]<-crossdot$x-0.3
            asd2_a$y_new[2]<-a2*asd2_a$gg[2]+b2
            asd2<-rbind(asd2,asd2_a)

            opt.inclusion<-opt.interval[as.integer(crossdot$x)]+(crossdot$x-as.integer(crossdot$x))*(opt.interval[as.integer(crossdot$x)+1]-opt.interval[as.integer(crossdot$x)])

            la <- list(a = as.numeric(format(a1, digits = 2)),
                       b = as.numeric(format(b1, digits = 2)),
                       dose = as.numeric(format(opt.inclusion, digits = 2)),
                       r2 = as.numeric(format(data_err_test1$R2, digits = 2)))
            #eq1<-substitute(italic(y) ==  a %.% (italic(x)-dose) + b~","~italic(R)^2~"="~r2,la)
            lb <- list(a = as.numeric(format(a2, digits = 2)),
                       b = as.numeric(format(b2, digits = 2)),
                       dose = as.numeric(format(opt.inclusion, digits = 2)),
                       r2 = as.numeric(format(data_err_test2$R2, digits = 2)))
            #eq2<-substitute(italic(y) ==  a %.% (italic(x)-dose) + b~","~italic(R)^2~"="~r2,lb)
            eq1<-""
            eq2<-""
            p<-p + geom_line(data = asd1,
                             aes(gg, y_new),
                             color = spline.color,
                             linewidth = spline.size/2,
                             linetype=spline.type,
                             lineend = "round")+
              geom_line(data = asd2,
                        aes(gg, y_new),
                        color = spline.color,
                        linewidth = spline.size/2,
                        linetype=spline.type,
                        lineend = "round")+
              annotate(geom="segment",x = crossdot$x, xend = crossdot$x,
                       size=vline.size,linetype=vline.type,
                       lineend = "round",
                       y =-Inf, yend = a1*crossdot$x+b1,
                       color=vline.color)+
              geom_point(aes(x=crossdot$x,y=a1*crossdot$x+b1),
                         shape=optimal.pt.type,
                         size=optimal.pt.size,color=optimal.pt.color)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=crossdot$x,y=a1*crossdot$x+b1),
                                       aes(x=crossdot$x,y=-Inf,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                       arrow = arrow(length = unit(0.04, "npc")),family="serif",size=2.5,
                                       box.padding = 1)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=(crossdot$x+min(asd1$gg))/2,
                                                       y=(a1*crossdot$x+b1+0)*0.6),
                                       family="serif",size=2,arrow = arrow(length = unit(0.03, "npc")),
                                       aes(x = opt.dose, y = y, label = as.character(as.expression(eq1))),parse = TRUE)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=(crossdot$x+max(asd2$gg))/2,
                                                       y=(a2*crossdot$x+b2+0)*0.4),
                                       family="serif",size=2,arrow = arrow(length = unit(0.03, "npc")),
                                       aes(x = opt.dose, y = y, label = as.character(as.expression(eq2))),parse = TRUE)
          }
        }
      }
    }
    }
  }
  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-microchatParamStatobj$data_err$group%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p<-p+ scale_x_discrete(labels = xlabname)
  }

  p<-p+theme(aspect.ratio = 1,
             title =  element_text(size = 7,face="bold.italic"),
             axis.title = element_text(size = 6,face="bold"),
             axis.text = element_text(size = 4,face="bold"))

  if (method %in% c("anova","kruskal.test")) {
    y.offset<-1/20*max(data_poi$value)
  }

  if (method %in% c("t.test","wilcox.test")) {
    y.offset<-max(data.seg$y3)-max(data_poi$value)+y.ratio*2/3
  }

  if (panel.border=="all") p<-p+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = -Inf+1, xend = -Inf+1)+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = Inf, yend = Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = Inf+1, xend = Inf+1)

  if (panel.border=="axis") p<-p+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = -Inf+1, xend = -Inf+1)

  if (panel.border=="normal") p<-p+theme(
    panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5)
  )
  if (panel.border=="none") p<-p

  if (axis.x.angle.adjust) p<-p+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

  if (!is.null(mytheme)) p<-p+mytheme
  ggsave(paste(export_path,"/Parameter (",index,")_boxplot.pdf",sep = ""),
         units = "cm",
         width = 21/3,
         height = 21*p$theme$aspect.ratio/3,
         p)
  cat("Parametric properities boxplot has been exported. Please check it.","\n")
  cat("----------------------------------------------------------------------","\n")
  return(p)
}

"calcorderx"<-function (data.alpha, y.ratio, data_poi)
{
  data.seg <- data.alpha
  data.seg <- subset(data.seg, select = c(index, group1, group2,
                                          p.value, sig, yvalue1, yvalue2, order, error))

  data.interset<-list()
  for (tt in 1:nrow(data.seg)) {
    data.interset[[tt]]<-c(match(data.seg$group1[tt],unique(data_poi$group)):match(data.seg$group2[tt],unique(data_poi$group)))
  }

  data.seg$yvalue1[1]<-data.seg$yvalue1[1]+y.ratio
  data.seg$yvalue2[1]<-data.seg$yvalue2[1]+y.ratio

  if (nrow(data.seg)>=2) {
    for (k in 2:nrow(data.seg)) {
      if ( (data.seg$group1[k] %in% data.seg$group1[1:(k-1)]) | (data.seg$group1[k] %in% data.seg$group2[1:(k-1)])) {
        rownum1<- match(data.seg$group1[k], data.seg$group1[1:(k-1)])
        rownum2<- match(data.seg$group1[k], data.seg$group2[1:(k-1)])

        if (!is.na(rownum1)) {
          if (!is.na(rownum2)) {
            if (rownum1>rownum2) data.seg$yvalue1[k]<-data.seg$yvalue1[rownum1]+y.ratio * 2
            if (rownum1<rownum2) data.seg$yvalue1[k]<-data.seg$yvalue2[rownum2]+y.ratio * 2
          } else {
            data.seg$yvalue1[k]<-data.seg$yvalue1[rownum1]+y.ratio * 2
          }
        } else {
          if (!is.na(rownum2)) {
            data.seg$yvalue1[k]<-data.seg$yvalue2[rownum2]+y.ratio * 2
          } else {
            data.seg$yvalue1[k]<-data.seg$yvalue1[k]
          }
        }

      }

      if ((data.seg$group2[k] %in% data.seg$group1[1:(k-1)]) | (data.seg$group2[k] %in% data.seg$group2[1:(k-1)])) {
        rownum3<- match(data.seg$group2[k], data.seg$group1[1:(k-1)])
        rownum4<- match(data.seg$group2[k], data.seg$group2[1:(k-1)])

        if (!is.na(rownum3)) {
          if (!is.na(rownum4)) {
            if (rownum3>rownum4) data.seg$yvalue2[k]<-data.seg$yvalue1[rownum3]+y.ratio * 2
            if (rownum3<rownum4) data.seg$yvalue2[k]<-data.seg$yvalue2[rownum4]+y.ratio * 2
          } else {
            data.seg$yvalue2[k]<-data.seg$yvalue1[rownum3]+y.ratio * 2
          }
        } else {
          if (!is.na(rownum4)) {
            data.seg$yvalue2[k]<-data.seg$yvalue2[rownum4]+y.ratio * 2
          } else {
            data.seg$yvalue2[k]<-data.seg$yvalue2[k]
          }
        }

      }

    }

    data.seg$order1 <- match(data.seg$group1, unique(data_poi$group))
    data.seg$order2 <- match(data.seg$group2, unique(data_poi$group))
    data.seg$y1 <- data.seg$yvalue1
    data.seg$y2 <- data.seg$yvalue2
    data.seg$y3 <- data.seg$y1 + y.ratio
    data.seg$y4 <- data.seg$y2 + y.ratio

    for (kk in 1:nrow(data.seg)) {
      if (data.seg$y3[kk] < data.seg$y4[kk]) {
        data.seg$y3[kk] <- data.seg$y4[kk]
      }
      else {
        data.seg$y3[kk] <- data.seg$y3[kk]
        data.seg$y4[kk] <- data.seg$y3[kk]
      }
    }

    data.intersetx<-list()
    for (tt in 1:nrow(data.seg)) {
      data.intersetx[[tt]]<-c(data.seg$order1[tt],data.seg$order2[tt])
      names(data.intersetx[[tt]])<-c(data.seg$group1[tt],data.seg$group2[tt])
    }

    ###判断左边还是右边会与其他comparison相交
    for (k in 2:nrow(data.seg)) {
      for (t in 1:(k-1)) {
        matchlen<-intersect(data.interset[[k]], data.interset[[t]])%>%length()

        if (matchlen>=2) {
          data.seg$y3[k]<- data.seg$y3[k-1]+y.ratio * 2
          data.seg$y4[k]<- data.seg$y4[k-1]+y.ratio * 2
        }

        order.n<-match(intersect(data.intersetx[[k]], data.interset[[t]]),data.intersetx[[k]])
        if (length(order.n)>0) {
          order.ns<-data.intersetx[[k]][order.n]%>%names()
          if (data.seg$group1[k]==order.ns) data.seg$y1[k]<-data.seg$y3[t]+y.ratio
          if (data.seg$group2[k]==order.ns) data.seg$y2[k]<-data.seg$y3[t]+y.ratio
        }

      }
    }
  } else {
    data.seg$order1 <- match(data.seg$group1, unique(data_poi$group))
    data.seg$order2 <- match(data.seg$group2, unique(data_poi$group))
    data.seg$y1 <- data.seg$yvalue1
    data.seg$y2 <- data.seg$yvalue2
    data.seg$y3 <- data.seg$y1 + y.ratio
    data.seg$y4 <- data.seg$y2 + y.ratio
  }


  return(data.seg)
}
"plotParametr"<-function(params,
                         color_taxa=NULL,
                         color_group=NULL,
                         bubble.color.taxa=TRUE,
                         bubble.max.size=9,
                         bubble.shape=19,
                         bubble.alpha=1,
                         aspect.ratio=0.8,
                         color_background=NULL,
                         color_border=NULL,
                         export_path="microbial composition") {
  export_path<-paste0(export_path,"/bubble plot")

  dir.create(export_path, recursive = TRUE)

  taxa_num<-ncol(params)
  genus_top<-t(params)%>%data.frame()

  if (bubble.color.taxa) {
    cat("\n","color taxa according to the params 'color_taxa'")
    genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
    genus_topxx<-reshape2::melt(genus_topxx)

    color_taxa<-c(color_taxa[1:taxa_num])
    names(color_taxa)<-unique(genus_topxx$tax)

    p1<-ggplot(data=genus_topxx, mapping=aes(x=variable,y=tax,color=tax))+
      geom_point(stat= "identity",aes(size=value),
                 shape=bubble.shape,
                 alpha=bubble.alpha,show.legend = TRUE)+
      scale_size(range = c(0.1, bubble.max.size),guide=FALSE)+
      scale_color_manual(values=color_taxa)+
      theme(aspect.ratio = aspect.ratio,
            axis.ticks.length = unit(0.2,"lines"),
            axis.ticks = element_line(color='black'),
            panel.background = element_rect(fill = "grey90"),
            #axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
            axis.text.y=element_text(colour='black',size=10,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 10,
                                     angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))
  } else {
    cat("\n","color taxa according to the params 'color_group'")

    genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
    genus_topxx<-reshape2::melt(genus_topxx)

    samplename<-unique(genus_topxx$variable)
    groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

    color_group<-color_group[1:length(groupname)]
    names(color_group)<-groupname

    color_group1<-color_group%>%data.frame()
    color_sample<-randomcolor(length(samplename))

    ###color samples according to groups
    for (ttks in 1:nrow(color_group1)) {
      color_sample[which(substr(samplename,start = 1, stop = 2)==rownames(color_group1)[ttks])]<-color_group1[ttks,]
    }
    names(color_sample)<-samplename

    p1<-ggplot(data=genus_topxx, mapping=aes(x=variable,y=tax,color=factor(variable)))+
      geom_point(stat= "identity",aes(size=value),
                 shape=bubble.shape,
                 alpha=bubble.alpha,show.legend = TRUE)+
      scale_size(range = c(0.1, bubble.max.size),guide=FALSE)+
      scale_color_manual(values=color_sample)+
      theme(aspect.ratio = aspect.ratio,
            axis.ticks.length = unit(0.2,"lines"),
            axis.ticks = element_line(color='black'),
            panel.background = element_rect(fill = "grey90"),
            #axis.line = element_line(colour = "black"),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(colour='black',size=10,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 10,
                                     angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))
  }

  if (!is.null(color_background)) {
    p1<-p1+
      theme(panel.background=element_rect(fill = color_background))
  }
  if (is.null(color_border)) {
    p1<-p1
  } else {
    p1<-p1+theme(panel.border = element_rect(fill=NA,colour = color_border))
  }
  p1<-p1+theme(aspect.ratio = ncol(params)/nrow(params))
  ggsave(paste(export_path,"/","bubble plot",".pdf",sep = ""),
         width = nrow(params)/2,
         height = ncol(params)/2,
         p1)
  return(p1)
}

"plotMicrochatParamHeatmap"<-function(microchatParamobj,
                                    rescale=FALSE,
                                    standardlization=c("0-1","scale","center","normal"),
                                    ###0-1 standardlization
                                    sample_type=c("allsamples","groups"),
                                    color_taxa=NULL,
                                    color_group=NULL,
                                    color.heatmap=NULL,  ###two colors at least
                                    heatmap.add.line=TRUE,
                                    heatmap.top.class=5,
                                    heatmap.sample.angle=90,
                                    heatmap.taxa.angle=0,
                                    export_path="microbial composition") {
  sample_type=match.arg(sample_type)
  standardlization=match.arg(standardlization)
  dir.create(export_path, recursive = TRUE)


  if (heatmap.top.class>(ncol(microchatParamobj$param_table)-2)) heatmap.top.class<-as.integer((ncol(microchatParamobj$param_table)-2)/2) else heatmap.top.class
  if (!rescale) {
    if (sample_type=="allsamples") {
      genus_top<-microchatParamobj$param_table
      genus_top<-genus_top[,2:ncol(genus_top)]
      library(tibble)
      genus_top<-tibble::column_to_rownames(genus_top,var="sample")
      genus_top<-t(genus_top)%>%data.frame()
    }
    if (sample_type=="groups") {
      genus_top<-microchatParamobj$param_table
      genus_top<-genus_top[,c(1,3:ncol(genus_top))]
      genus_top<-genus_top %>%
        group_by(group) %>%
        summarise_all(sum)
      genus_top<-tibble::column_to_rownames(genus_top,var="group")
      genus_top<-t(genus_top)%>%data.frame()
    }
  } else {
    if (sample_type=="allsamples") {
      if (standardlization=="0-1") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,2:ncol(genus_top)]
        library(tibble)
        genus_top<-tibble::column_to_rownames(genus_top,var="sample")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<- scale1dt(genus_top)
      }
      if (standardlization=="normal") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,2:ncol(genus_top)]
        library(tibble)
        genus_top<-tibble::column_to_rownames(genus_top,var="sample")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<-scale4dt(genus_top)
      }
      if (standardlization=="scale") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,2:ncol(genus_top)]
        library(tibble)
        genus_top<-tibble::column_to_rownames(genus_top,var="sample")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<- scale2dt(genus_top)
      }
      if (standardlization=="center") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,2:ncol(genus_top)]
        library(tibble)
        genus_top<-tibble::column_to_rownames(genus_top,var="sample")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<- scale3dt(genus_top)
      }
    }
    if (sample_type=="groups") {

      if (standardlization=="0-1") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,c(1,3:ncol(genus_top))]
        genus_top<-genus_top %>%
          group_by(group) %>%
          summarise_all(sum)
        genus_top<-tibble::column_to_rownames(genus_top,var="group")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<-scale1dt(genus_top)
      }

      if (standardlization=="normal") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,c(1,3:ncol(genus_top))]
        genus_top<-genus_top %>%
          group_by(group) %>%
          summarise_all(sum)
        genus_top<-tibble::column_to_rownames(genus_top,var="group")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<-scale4dt(genus_top)
      }


      if (standardlization=="scale") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,c(1,3:ncol(genus_top))]
        genus_top<-genus_top %>%
          group_by(group) %>%
          summarise_all(sum)
        genus_top<-tibble::column_to_rownames(genus_top,var="group")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<-scale2dt(genus_top)
      }
      if (standardlization=="center") {
        genus_top<-microchatParamobj$param_table
        genus_top<-genus_top[,c(1,3:ncol(genus_top))]
        genus_top<-genus_top %>%
          group_by(group) %>%
          summarise_all(sum)
        genus_top<-tibble::column_to_rownames(genus_top,var="group")
        genus_top<-t(genus_top)%>%data.frame()
        genus_top<-scale3dt(genus_top)
      }
    }
  }

  net<-genus_top%>%as.matrix()

  samplename<-colnames(net)
  groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

  ###color groups
  if (is.null(color_group)) {
    color_group<-randomcolor(length(groupname))
  } else {
    if (length(color_group)<length(groupname)) {
      stop("There are no enough colors to fill all groups !!!")
    } else {
      color_group<-color_group[1:length(groupname)]
    }
  }

  names(color_group)<-groupname

  color_group1<-color_group%>%data.frame()
  color_sample<-randomcolor(length(samplename))

  ###color samples according to groups
  for (ttks in 1:nrow(color_group1)) {
    color_sample[which(substr(samplename,start = 1, stop = 2)==rownames(color_group1)[ttks])]<-color_group1[ttks,]
  }
  names(color_sample)<-samplename

  ###color taxa
  taxa_name<-rownames(net)
  if (is.null(color_taxa)) {
    color_taxa<-randomcolor(length(taxa_name))
  } else {
    if (length(color_taxa)<length(taxa_name)) {
      stop("There are no enough colors to fill all taxa !!!")
    } else {
      color_taxa<-color_taxa[1:length(taxa_name)]
    }
  }
  names(color_taxa)<-taxa_name

  ###color heatmap blocks
  if (is.null(color.heatmap)) {
    color.heatmap.use = grDevices::colorRampPalette((
      RColorBrewer::brewer.pal(n = 9,name = "GnBu")
    ))(100)
  } else {
    if (length(color.heatmap)<2) {
      stop("There are no enough colors to fill heatmap blocks !!!")
    } else {
      color.heatmap.use = grDevices::colorRampPalette((
        color.heatmap
      ))(100)
    }
  }

  df1 <- data.frame(group = colnames(net))
  rownames(df1) <- colnames(net)
  color.use1 = color_sample

  col_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df1, col = list(group = color.use1),
                                                      which = "column", show_legend = TRUE, show_annotation_name = FALSE,
                                                      simple_anno_size = grid::unit(0.2, "cm"),
                                                      annotation_legend_param=list(title_gp = grid::gpar(fontsize = 8,
                                                                                                         fontfamily="serif",
                                                                                                         fontface = "plain"),
                                                                                   title_position = "leftcenter-rot",
                                                                                   border = NA,
                                                                                   legend_height = unit(20, "mm"),
                                                                                   labels_gp = grid::gpar(fontsize = 8,fontfamily="serif"),
                                                                                   grid_width = unit(2, "mm")))

  df2 <- data.frame(group = rownames(net))
  rownames(df2) <- rownames(net)
  color.use2 = color_taxa

  row_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df2, col = list(group = color.use2),
                                                      which = "row", show_legend = TRUE, show_annotation_name = FALSE,
                                                      simple_anno_size = grid::unit(0.2, "cm"),
                                                      annotation_legend_param=list(title_gp = grid::gpar(fontsize = 8,
                                                                                                         fontfamily="serif",
                                                                                                         fontface = "plain"),
                                                                                   title_position = "leftcenter-rot",
                                                                                   border = NA,
                                                                                   legend_height = unit(20, "mm"),
                                                                                   labels_gp = grid::gpar(fontsize = 8,fontfamily="serif"),
                                                                                   grid_width = unit(2, "mm")))

  if (heatmap.add.line) {
    ha1 = ComplexHeatmap::rowAnnotation(
      Strength = ComplexHeatmap::anno_barplot(rowSums(abs(net)),
                                              border = FALSE, gp = grid::gpar(
                                                fill = color.use2, col = color.use2)),
      foo = ComplexHeatmap::anno_lines(rowSums(abs(net)),
                                       gp = grid::gpar(col = "red"),
                                       add_points = TRUE,smooth = TRUE, axis = FALSE,
                                       extend = 0.2,border = FALSE,
                                       pt_gp = grid::gpar(fill = color.use2,col = color.use2),
                                       height = unit(1, "cm")),
      show_annotation_name = FALSE)
    ha2 = ComplexHeatmap::HeatmapAnnotation(
      #pt = anno_points(colSums(abs(net[1:heatmap.top.class,])), gp = gpar(fill = color.use1,col = color.use1),height = unit(1, "cm")),
      foo = ComplexHeatmap::anno_lines(colSums(abs(net[1:heatmap.top.class,])),
                                       gp = grid::gpar(col = "red"),
                                       add_points = TRUE,smooth = TRUE, axis = FALSE,
                                       extend = 0.2,border = FALSE,
                                       pt_gp = grid::gpar(fill = color.use1,col = color.use1),
                                       height = unit(1, "cm")),
      Strength = ComplexHeatmap::anno_barplot(colSums(abs(net)),
                                              border = FALSE, gp = grid::gpar(
                                                fill = color.use1, col = color.use1)),
      show_annotation_name = FALSE)
  } else {
    ha1 = ComplexHeatmap::rowAnnotation(
      Strength = ComplexHeatmap::anno_barplot(rowSums(abs(net)),
                                              border = FALSE, gp = grid::gpar(
                                                fill = color.use2, col = color.use2)),
      show_annotation_name = FALSE)
    ha2 = ComplexHeatmap::HeatmapAnnotation(
      Strength = ComplexHeatmap::anno_barplot(colSums(abs(net)),
                                              border = FALSE, gp = grid::gpar(
                                                fill = color.use1, col = color.use1)),
      show_annotation_name = FALSE)

  }


  p1<-ComplexHeatmap::Heatmap(net,
                              col = color.heatmap.use,
                              na_col = "white",
                              name = paste("Value",sep = ""),
                              bottom_annotation = col_annotation,
                              left_annotation = row_annotation,
                              top_annotation = ha2,
                              right_annotation = ha1,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_names_side = "left",
                              row_names_rot = heatmap.taxa.angle,
                              row_names_gp = grid::gpar(fontsize = 6,
                                                        fontfamily="serif"),
                              column_names_gp = grid::gpar(fontsize = 6,
                                                           fontfamily="serif"),
                              column_names_rot = heatmap.sample.angle,
                              heatmap_legend_param = list(nrow=3,
                                                          title_gp = grid::gpar(fontsize = 8,
                                                                                fontfamily="serif",
                                                                                fontface = "plain"),
                                                          title_position = "leftcenter-rot",
                                                          border = NA,
                                                          legend_height = unit(20, "mm"),
                                                          labels_gp = grid::gpar(fontsize = 8,fontfamily="serif"),
                                                          grid_width = unit(2, "mm")))


  message("\n","The heatmap need to be saved manually.")
  return(p1)
}


"calcSegment3Order" <- function(data.alpha,y.ratio,data_poi) {
  data.er<-subset(data.alpha, select = c(index,group1,group2,max1,max2,sig,order,error))
  data.er$order1<-match(data.er$group1,unique(data_poi$group))
  data.er$order2<-match(data.er$group2,unique(data_poi$group))

  data.er$y1<-data.er$max1
  data.er$y2<-data.er$max2
  data.er$y3<-data.er$y1+y.ratio
  data.er$y4<-data.er$y2+y.ratio

  data.er$y1<-data.er$y1+y.ratio
  data.er$y2<-data.er$y2+y.ratio
  data.er$y3<-data.er$y3+y.ratio
  data.er$y4<-data.er$y4+y.ratio

  for (tt in 1:nrow(data.er)) {
    x1<-data.er$order1[tt]
    x2<-data.er$order2[tt]

    y1.use<-data.er$y1[tt]
    y2.use<-data.er$y2[tt]
    y3.use<-data.er$y3[tt]
    y4.use<-data.er$y4[tt]

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

    data.er$y1[tt]<-y1.use
    data.er$y2[tt]<-y2.use
    data.er$y3[tt]<-y3.use
    data.er$y4[tt]<-y4.use

  }

  order<-unique(data.er$order)
  for (index in order[2:length(order)]) {
    data.er[which(data.er$order==index),]$y1<-max(data.er[which(data.er$order==(index-1)),]$y3)+y.ratio
    data.er[which(data.er$order==index),]$y3<-data.er[which(data.er$order==index),]$y1+y.ratio
  }
  return(data.er)

}

"plotMicrochatParamMutiBoxplot" <- function(microchatParamobj,
                                            geom.line=TRUE,
                                            geom.line.size=0.5,
                                            ylim.fold=1.1,
                                            xlabname=NULL,
                                            yaxis.italic=TRUE,
                                            strictmod=TRUE,
                                            panel.border=c("none","normal","all","axis"),
                                            method="anova",
                                            seg=FALSE,
                                            comparison=my_comparisons,
                                            color_group=colorCustom(5,pal = "ywbu"),
                                            color_backgroud="grey90",
                                            axis.x.angle.adjust=FALSE,
                                            mytheme=NULL,
                                            ncol=3,
                                            spline.params=NULL,
                                            export_path="ss21/microbial parameteric analysis/liver_gene") {
  if (dev.cur() != 1) {
    dev.off()
  }
  select.index<-microchatParamobj$all.index
  paramfile.select<-microchatParamobj$paramfile.select
  export_path<-paste0(export_path,"/",paramfile.select)
  pp<-list()
  height.sel<-length(ncol_layout(length(select.index)))*3
  width.sel<-ncol_layout(length(select.index))[[1]]%>%length()*3
  #t=select.index[8]
  for (t in select.index) {
    microchatParamStatobj<-calcMicrochatParamStat(microchatParamobj,
                                                  select.index=t,
                                                  strictmod=strictmod,
                                                  method=method,
                                                  comparison=comparison,
                                                  export_path=export_path)

    p<-plotMicrochatParamBoxplot(microchatParamStatobj,
                                 geom.line=geom.line,
                                 geom.line.size=geom.line.size,
                                 ylim.fold=ylim.fold,
                                 xlabname=xlabname,
                                 yaxis.italic=yaxis.italic,
                                 errorbar.line.add=TRUE,
                                 errorbar.point.size=0.1,
                                 y.point.adj=0.1,
                                 seg=seg,
                                 panel.border=panel.border,
                                 color_group=color_group,
                                 color_backgroud=color_backgroud,
                                 axis.x.angle.adjust=axis.x.angle.adjust,
                                 mytheme=mytheme,
                                 spline.params=spline.params,
                                 export_path=export_path)
    pp[[t]]<-p
  }

  library(patchwork)
  pm<-wrap_plots(pp,ncol=ncol)

  rownum<-length(select.index)/ncol
  if (rownum>as.integer(rownum)) nrows<-as.integer(rownum)+1 else nrows<-as.integer(rownum)

  ggsave(paste(export_path,"/Muti-parameter","_boxplot.pdf",sep = ""),
         units = "cm",
         width = 7*ncol,height = 7*nrows,
         pm)
  message("Muti-parametric properities boxplot has been exported. Please check it.")

  return(pm)
}

"calcMicrochatParamTable" <- function(microchatParamobj,
                                      strictmod=TRUE,
                                      method="anova",
                                      comparison=my_comparisons,
                                      mean_dec=2,
                                      se_dec=2,
                                      export_path="microbial parameteric analysis") {
  paramfile.select<-microchatParamobj$paramfile.select
  export_path<-paste(export_path,"/",paramfile.select,sep = "")
  dir.create(export_path, recursive = TRUE)
  if (class(microchatParamobj)[1]!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  suppressMessages(library(reshape2))
  suppressMessages(library(microchat))

  data_err<-microchatParamobj$statistics
  data_err$summa<-paste(sprintf(paste0("%0.",mean_dec,"f"), data_err$mean),sprintf(paste0("%0.",se_dec,"f"), data_err$se),sep = "±")

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
        if (sigtab$sig[tk]=="**") sigtab$sig[tk]<-"##"
        if (sigtab$sig[tk]=="***") sigtab$sig[tk]<-"###"
        if (sigtab$sig[tk]=="****") sigtab$sig[tk]<-"####"
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
    anova_p<-data.frame()
    t=select.index[1]
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
        anova_p<-rbind(anova_p,microchatParamStatobj$param.stats)
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

    param_tab$ANOVA_p<- keep.decmi(anova_p[match(rownames(param_tab),anova_p$index),]$p.value)
    ###add opc prediction
    params_table<-microchatParamobj$param_table
    params_table<-params_table[,-1]
    params_table<-column_to_rownames(params_table,var = "sample")%>%data.frame()

    colname<-colnames(params_table)
    trt_id <-unique(colname)

    ttk<-sapply(trt_id,function(x){
      match(x,colnames(params_table)) #grep(x,colnames(params_table))
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
  class(param_tab) <- c("microchat","data.frame")
  return(param_tab)
}

"plotMicrochatParamComplexBoxplot" <- function(submchat,
                                               geom.line=TRUE,
                                               geom.line.size=0.5,
                                               ylim.fold=1.1,
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
                                  geom.line=geom.line,
                                  yaxis.italic=yaxis.italic,
                                  strictmod=strictmod,
                                  method=method,
                                  comparison=comparison,
                                  color_group=color_group,
                                  export_path=export_path2)
  }
}



"plotMicrochatParamBarplot" <- function(microchatParamStatobj,
                                        ylim.fold=1.1,
                                        xlabname=NULL,
                                        yaxis.italic=TRUE,
                                        errorbar.pos.adj=TRUE,
                                        errorbar.line.add=FALSE,
                                        seg=TRUE,
                                        add.spline=TRUE,
                                        panel.border=c("none","normal","all","axis"),
                                        errorbar.point.size=0.1,
                                        y.point.adj=NULL,
                                        color_group=colorCustom(5,pal = "gygn"),
                                        color_backgroud="grey90",
                                        axis.x.angle.adjust=FALSE,
                                        mytheme=NULL,
                                        spline.params=NULL,
                                        export_path="microbial diversity analysis") {
  panel.border<-match.arg(panel.border)
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

  if (method %in% c("anova","kruskal.test")) {
    opc.data <- microchatParamStatobj$opc.data
    o.data <- opc.data[which(opc.data$index == index), ]
    gnum <- ncol(o.data)
    gnum/2
    liner <- paste(colnames(o.data[gnum/2 + 2]), " effect: ",
                   o.data[gnum/2 + 2], sep = "")
    quad <- paste(colnames(o.data[gnum/2 +3]), " effect: ",
                  o.data[gnum/2 +  3], sep = "")
    subtitile <- paste(liner, "\n",quad, sep = "")
}

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

    if (!seg) data.seg <- calcSegmentOrder(data.alpha, y.ratio,
                                           data_poi)
    if (seg) data.seg <-calcorderx(data.alpha, y.ratio,
                                   data_poi)
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

      data.seg$y1[tt]<-y1.use
      data.seg$y2[tt]<-y2.use
      data.seg$y3[tt]<-y3.use
      data.seg$y4[tt]<-y4.use

    }

    for (tt in 1:nrow(data.seg)) {
      x1<-data.seg$order1[tt]
      x2<-data.seg$order2[tt]

      y1.use<-data.seg$y1[tt]
      y2.use<-data.seg$y2[tt]
      y3.use<-data.seg$y3[tt]
      y4.use<-data.seg$y4[tt]

      label.use=data.seg$sig[tt]

      p<-p+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x1, y = y1.use, ymin = y1.use, ymax = y1.use)+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x2, y = y2.use, ymin = y2.use, ymax = y2.use)+

        annotate("segment",color="grey50", x = x1, y = y1.use, xend = x1, yend = y3.use)+
        annotate("segment",color="grey50", x = x2, y = y2.use, xend = x2, yend = y3.use)+
        annotate("segment",color="grey50", x = x1, y = y3.use, xend = x2, yend = y3.use)+
        annotate("text",family="serif",x = (x1+x2)/2, y = y3.use+y.ratio*2/3,label=label.use)

      if (errorbar.line.add) p<-p+annotate("segment", color="grey50",x = x1-0.05, y = y1.use, xend = x1+0.05, yend = y1.use)+
        annotate("segment",color="grey50", x = x2-0.05, y = y2.use, xend = x2+0.05, yend = y2.use)

    }

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


    if (yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="bold.italic", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

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
                  size = 3,color = "black",family = "serif", fontface="bold")


      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")

      if (yaxis.italic) p<-p+labs(y=index,
                                  title = paste(
                                    paste(method,": p = ",
                                          keep.decmi(data.alpha$p.value),sep=""),"\n",
                                    subtitile,sep = ""
                                  ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   title = paste(
                                     paste(method,": p = ",
                                           keep.decmi(data.alpha$p.value),sep=""),"\n",
                                     subtitile,sep = ""
                                   ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

      p<-ggplot() +geom_errorbar(data = data_err,
                                 aes(x = group, y = mean, group = index,
                                     ymin = mean-se, ymax = mean+se),
                                 colour = "grey50",width=.25)+
        geom_bar(data = data_err,aes(x = group,y = mean,fill = group),
                 stat = "identity",position = "dodge",
                 width = 0.7,colour = "white") +

        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 3,color = "black",family = "serif", fontface="bold")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")

      if (yaxis.italic) p<-p+labs(y=index,
                                  title = paste(
                                    paste(method,": p = ",
                                          keep.decmi(data.alpha$p.value),sep=""),"\n",
                                    subtitile,sep = ""
                                  ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   title = paste(
                                     paste(method,": p = ",
                                           keep.decmi(data.alpha$p.value),sep=""),"\n",
                                     subtitile,sep = ""
                                   ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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

  ##add fitted curve for  quadratic model
  if (method %in% c("anova","kruskal.test") & !is.null(spline.params)){
    opc.thres=spline.params$opc.thres
    spline.size = spline.params$spline.size
    spline.color=spline.params$spline.color
    spline.shape =spline.params$spline.shape
    spline.type=spline.params$spline.type
    opt.interval=spline.params$opt.interval
    optimal.pt.size=spline.params$optimal.pt.size
    optimal.pt.color=spline.params$optimal.pt.color
    optimal.pt.type=spline.params$optimal.pt.type

    vline.color=spline.params$vline.color
    vline.size=spline.params$vline.size
    vline.type=spline.params$vline.type

    csal<-strsplit(o.data$Linear,"( )")[[1]][2]
    lin.num<-substr(csal,start = 2,stop =nchar(csal)-1 )%>%as.numeric()

    csa<-strsplit(o.data$Quadratic,"( )")[[1]][2]
    qua.num<-substr(csa,start = 2,stop =nchar(csa)-1 )%>%as.numeric()
    if (!is.na(lin.num)) {
    if ("Cubic" %in% colnames(o.data)) {
      csac<-strsplit(o.data$Cubic,"( )")[[1]][2]
      cub.num<-substr(csac,start = 2,stop =nchar(csac)-1 )%>%as.numeric()

      all.opc.r<-c(lin.num,qua.num,cub.num)

    } else {
      all.opc.r<-c(lin.num,qua.num)
    }

    #拟合三次函数
      if ("Cubic" %in% colnames(o.data)) if (cub.num>opc.thres & max(all.opc.r)==cub.num) {
      data_err_test<-cub_equa_fitTest(microchatParamStatobj$data_err)

      asd<-data_err_test$asd
      opt.dose<-data_err_test$opt.dose
      a<-data_err_test$a
      b<-data_err_test$b
      c<-data_err_test$c
      d<-data_err_test$d
      opt.inclusion<-opt.interval[as.integer(opt.dose)]+(opt.dose-as.integer(opt.dose))*(opt.interval[as.integer(opt.dose)+1]-opt.interval[as.integer(opt.dose)])

      l <- list(d = as.numeric(format(d, digits = 2)),
                a = as.numeric(format(a, digits = 2)),
                b = as.numeric(format(b, digits = 2)),
                c = as.numeric(format(c, digits = 2)),
                dose = as.numeric(format(opt.inclusion, digits = 2)),
                r2 = as.numeric(format(data_err_test$R2, digits = 2)))
      #eq<-substitute(italic(y) ==  d %.% (italic(x)-dose)^3 + a %.% (italic(x)-dose)^2 + b %.% (italic(x)-dose) + c~","~italic(R)^2~"="~r2,l)
      eq<-""

      p<-p + geom_line(data = data.frame(x = data_err_test$x_pred,
                                         y = data_err_test$y_pred),
                       aes(x, y),
                       color = spline.color,
                       linewidth = spline.size/2,
                       linetype=spline.type,
                       lineend = "round")+
        annotate(geom="segment",x = opt.dose, xend = opt.dose, size=vline.size,linetype=vline.type,
                 lineend = "round",
                 y =-Inf, yend = d*opt.dose^3+ a*opt.dose^2+b*opt.dose+c,color=vline.color)+
        geom_point(aes(x=opt.dose,y=d*opt.dose^3+a*opt.dose^2+b*opt.dose+c),
                   shape=optimal.pt.type,
                   size=optimal.pt.size,color=optimal.pt.color)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(d*opt.dose^3+a*opt.dose^2+b*opt.dose+c)),
                                 aes(x=opt.dose,y=y,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                 arrow = arrow(length = unit(0.01, "npc")),family="serif",size=2.5,
                                 box.padding = 1)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)*0.8),
                                 family="serif",size=2,
                                 aes(x = opt.dose, y = -Inf, label = as.character(as.expression(eq))),parse = TRUE)
    }
    #拟合二次函数
    if (qua.num>opc.thres & max(all.opc.r)==qua.num) {
      data_err_test<-qua_equa_fitTest(microchatParamStatobj$data_err)
      asd<-data_err_test$asd
      opt.dose<-data_err_test$opt.dose
      a<-data_err_test$a
      b<-data_err_test$b
      c<-data_err_test$c
      opt.inclusion<-opt.interval[as.integer(opt.dose)]+(opt.dose-as.integer(opt.dose))*(opt.interval[as.integer(opt.dose)+1]-opt.interval[as.integer(opt.dose)])

      l <- list(a = as.numeric(format(a, digits = 3)),
                b = as.numeric(format(b, digits = 3)),
                c=as.numeric(format(c, digits = 3)),
                dose = as.numeric(format(opt.inclusion, digits = 2)),
                r2 = as.numeric(format(data_err_test$R2, digits = 3)))
      #eq<-substitute(italic(y) == a %.% (italic(x)-dose)^2 + b %.% (italic(x)-dose) + c~","~italic(R)^2~"="~r2,l)
      eq<-""

      p<-p + ggalt::geom_xspline(data = asd,
                                 linetype=spline.type,lineend = "round",
                                 size = spline.size, color = spline.color,
                                 aes(x = gg, y = y_new), spline_shape = spline.shape) +
        annotate(geom="segment",x = opt.dose, xend = opt.dose, size=vline.size,linetype=vline.type,
                 lineend = "round",
                 y =-Inf, yend =  a*opt.dose^2+b*opt.dose+c,color=vline.color)+
        geom_point(aes(x=opt.dose,y=a*opt.dose^2+b*opt.dose+c),
                   shape=optimal.pt.type,
                   size=optimal.pt.size,color=optimal.pt.color)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)),
                                 aes(x=opt.dose,y=y,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                 arrow = arrow(length = unit(0.01, "npc")),family="serif",size=2.5,
                                 box.padding = 1)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)*0.8),
                                 family="serif",size=2,
                                 aes(x = opt.dose, y = -Inf, label = as.character(as.expression(eq))),parse = TRUE)
    }
    #拟合线性模型
    if (lin.num>opc.thres & max(all.opc.r)==lin.num) {
      data_err_test<-linear_equa_fitTest(microchatParamStatobj$data_err)

      if (data_err_test$a>=0) {
        message("The fitted curve was drawn as user required. However, the slope of the fitted curve is higher than 0, so the fitted curve will not display.")
      } else {
        message("The fitted curve was drawn as user required. The slope of the fitted curve is lower than 0, so the fitted curve will be divided into two segments.")
        message("Namely, fitted curve for broken-line model is showing  !!!")

        slope<-c()
        for (t in 1:(nrow(data_err_test$asd)-1)) {
          slope<-c(slope,obtainSlope(data<-data_err_test$asd[t:(t+1),]))
        }
        slope_a<-slope[which(slope<0)]
        slope_a<-slope_a[which(slope_a==min(slope_a))]
        spot_start<-which(slope==slope_a)
        if (which(data_err_test$asd$mean==max(data_err_test$asd$mean))!=as.numeric(spot_start)) {
          spot_start<-which(data_err_test$asd$mean==max(data_err_test$asd$mean))
        } else {
          spot_start<-spot_start
        }
        if (spot_start ==2 | spot_start==3) {
          spot_end<-spot_start+1

          data_err_test1<-linear_equa_fitTest(microchatParamStatobj$data_err[1:spot_start,])
          data_err_test2<-linear_equa_fitTest(microchatParamStatobj$data_err[spot_end:nrow(microchatParamStatobj$data_err),],spot_end)

          crossdot<-getCrossCoord(data_err_test1,data_err_test2)

          if (crossdot$x>1 & crossdot$x <nrow(data_err_test$asd)) {
            asd1<-data_err_test1$asd
            a1<-data_err_test1$a
            b1<-data_err_test1$b
            asd2<-data_err_test2$asd
            asd2$gg<-(max(asd1$gg)+1):nrow(data_err_test$asd)
            a2<-data_err_test2$a
            b2<-data_err_test2$b

            asd1_a<-asd1[1:2,]
            asd1_a$gg[1]<-min(asd1$gg)-0.3
            asd1_a$y_new[1]<-a1*asd1_a$gg[1]+b1
            asd1_a$gg[2]<-crossdot$x+0.3
            asd1_a$y_new[2]<-a1*asd1_a$gg[2]+b1
            asd1<-rbind(asd1,asd1_a)

            asd2_a<-asd2[1:2,]
            asd2_a$gg[1]<-max(asd2$gg)+0.3
            asd2_a$y_new[1]<-a2*asd2_a$gg[1]+b2
            asd2_a$gg[2]<-crossdot$x-0.3
            asd2_a$y_new[2]<-a2*asd2_a$gg[2]+b2
            asd2<-rbind(asd2,asd2_a)

            opt.inclusion<-opt.interval[as.integer(crossdot$x)]+(crossdot$x-as.integer(crossdot$x))*(opt.interval[as.integer(crossdot$x)+1]-opt.interval[as.integer(crossdot$x)])

            la <- list(a = as.numeric(format(a1, digits = 2)),
                       b = as.numeric(format(b1, digits = 2)),
                       dose = as.numeric(format(opt.inclusion, digits = 2)),
                       r2 = as.numeric(format(data_err_test1$R2, digits = 2)))
            #eq1<-substitute(italic(y) ==  a %.% (italic(x)-dose) + b~","~italic(R)^2~"="~r2,la)
            lb <- list(a = as.numeric(format(a2, digits = 2)),
                       b = as.numeric(format(b2, digits = 2)),
                       dose = as.numeric(format(opt.inclusion, digits = 2)),
                       r2 = as.numeric(format(data_err_test2$R2, digits = 2)))
            #eq2<-substitute(italic(y) ==  a %.% (italic(x)-dose) + b~","~italic(R)^2~"="~r2,lb)
            eq1<-""
            eq2<-''
            p<-p + geom_line(data = asd1,
                             aes(gg, y_new),
                             color = spline.color,
                             linewidth = spline.size/2,
                             linetype=spline.type,
                             lineend = "round")+
              geom_line(data = asd2,
                        aes(gg, y_new),
                        color = spline.color,
                        linewidth = spline.size/2,
                        linetype=spline.type,
                        lineend = "round")+
              annotate(geom="segment",x = crossdot$x, xend = crossdot$x,
                       size=vline.size,linetype=vline.type,
                       lineend = "round",
                       y =-Inf, yend = a1*crossdot$x+b1,
                       color=vline.color)+
              geom_point(aes(x=crossdot$x,y=a1*crossdot$x+b1),
                         shape=optimal.pt.type,
                         size=optimal.pt.size,color=optimal.pt.color)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=crossdot$x,y=a1*crossdot$x+b1),
                                       aes(x=crossdot$x,y=-Inf,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                       arrow = arrow(length = unit(0.04, "npc")),family="serif",size=2.5,
                                       box.padding = 1)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=(crossdot$x+min(asd1$gg))/2,
                                                       y=(a1*crossdot$x+b1+0)*0.6),
                                       family="serif",size=2,arrow = arrow(length = unit(0.03, "npc")),
                                       aes(x = opt.dose, y = y, label = as.character(as.expression(eq1))),parse = TRUE)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=(crossdot$x+max(asd2$gg))/2,
                                                       y=(a2*crossdot$x+b2+0)*0.4),
                                       family="serif",size=2,arrow = arrow(length = unit(0.03, "npc")),
                                       aes(x = opt.dose, y = y, label = as.character(as.expression(eq2))),parse = TRUE)
          }
        }
      }
    }
    }
  }

  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-microchatParamStatobj$data_err$group%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p<-p+ scale_x_discrete(labels = xlabname)
  }

  p<-p+theme(aspect.ratio = 1,
             title =  element_text(size = 7,face="bold.italic"),
             axis.title = element_text(size = 6,face="bold"),
             axis.text = element_text(size = 4,face="bold"))

  if (method %in% c("anova","kruskal.test")) {
    y.offset<-1/20*max(data_poi$value)
  }

  if (method %in% c("t.test","wilcox.test")) {
    y.offset<-max(data.seg$y3)-max(data_poi$value)+y.ratio*2/3
  }

  if (panel.border=="all") p<-p+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = -Inf+1, xend = -Inf+1)+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = Inf, yend = Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = Inf+1, xend = Inf+1)

  if (panel.border=="axis") p<-p+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = -Inf+1, xend = -Inf+1)

  if (panel.border=="normal") p<-p+theme(
    aspect.ratio = 1,
    title =  element_text(size = 7,face="bold.italic"),
    axis.title = element_text(size = 6,face="bold"),
    axis.text = element_text(size = 4,face="bold"),
    panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5)
  )
  if (panel.border=="none") p<-p

  if (axis.x.angle.adjust) p<-p+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
  if (!is.null(mytheme)) p<-p+mytheme
  ggsave(paste(export_path,"/Parameter (",index,")_barplot.pdf",sep = ""),
         units = "cm",
         width = 21/3,
         height = 21*p$theme$aspect.ratio/3,
         p)
  cat("Parametric properities barplot has been exported. Please check it.","\n")
  cat("----------------------------------------------------------------------","\n")
  return(p)

}

"plotMicrochatParamLineplot" <- function(microchatParamStatobj,
                                         ylim.fold=1.1,
                                         xlabname=NULL,
                                         yaxis.italic=TRUE,
                                         errorbar.pos.adj=TRUE,
                                         errorbar.line.add=FALSE,
                                         seg=TRUE,
                                         panel.border=TRUE,
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
  if (method %in% c("anova","kruskal.test")) {
    opc.data <- microchatParamStatobj$opc.data
    o.data <- opc.data[which(opc.data$index == index), ]
    gnum <- ncol(o.data)
    gnum/2
    liner <- paste(colnames(o.data[gnum/2 + 2]), " effect: ",
                   o.data[gnum/2 + 2], sep = "")
    quad <- paste(colnames(o.data[gnum/2 +3]), " effect: ",
                  o.data[gnum/2 +  3], sep = "")
    subtitile <- paste(liner, quad, sep = " ")
  }

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


    if (!seg) data.seg <- calcSegmentOrder(data.alpha, y.ratio,
                                           data_poi)
    if (seg) data.seg <-calcorderx(data.alpha, y.ratio,
                                   data_poi)
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

      data.seg$y1[tt]<-y1.use
      data.seg$y2[tt]<-y2.use
      data.seg$y3[tt]<-y3.use
      data.seg$y4[tt]<-y4.use

    }

    for (tt in 1:nrow(data.seg)) {
      x1<-data.seg$order1[tt]
      x2<-data.seg$order2[tt]

      y1.use<-data.seg$y1[tt]
      y2.use<-data.seg$y2[tt]
      y3.use<-data.seg$y3[tt]
      y4.use<-data.seg$y4[tt]

      label.use=data.seg$sig[tt]
      p<-p+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x1, y = y1.use, ymin = y1.use, ymax = y1.use)+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x2, y = y2.use, ymin = y2.use, ymax = y2.use)+

        annotate("segment",color="grey50", x = x1, y = y1.use, xend = x1, yend = y3.use)+
        annotate("segment",color="grey50", x = x2, y = y2.use, xend = x2, yend = y3.use)+
        annotate("segment",color="grey50", x = x1, y = y3.use, xend = x2, yend = y3.use)+
        annotate("text",family="serif",x = (x1+x2)/2, y = y3.use+y.ratio*2/3,label=label.use)

      if (errorbar.line.add) p<-p+annotate("segment", color="grey50",x = x1-0.05, y = y1.use, xend = x1+0.05, yend = y1.use)+
        annotate("segment",color="grey50", x = x2-0.05, y = y2.use, xend = x2+0.05, yend = y2.use)

    }

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


    if (yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="bold.italic", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      sig_label_new<-sig_label_new[order(sig_label_new$group),]
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

      p<-ggplot()+
        geom_line(data= data_err,size=0.5,colour = "grey50",aes(x=group, y=mean, group=index)) +
        geom_errorbar(data = data_err,
                      aes(x = group, y = mean, group = index,
                          ymin = mean-se, ymax = mean+se),
                      colour = "grey50",width=.25)+
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 5,color = "black",family = "serif", fontface="bold")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")


      if (yaxis.italic) p<-p+labs(y=index,
                                  subtitle=subtitile,
                                  title = paste(method,": p = ",
                                                keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   subtitle=subtitile,
                                   title = paste(method,": p = ",
                                                 keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      sig_label_new<-sig_label_new[order(sig_label_new$group),]
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

      p<-ggplot() +geom_line(data= data_err,size=0.5,colour = "grey50",
                             aes(x=group, y=mean, group=index))+
        geom_errorbar(data = data_err,
                      aes(x = group, y = mean, group = index,
                          ymin = mean-se, ymax = mean+se),
                      colour = "grey50",width=.25)+
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 5,color = "black",family = "serif", fontface="bold")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")
      if (yaxis.italic) p<-p+labs(y=index,
                                  subtitle=subtitile,
                                  title = paste(method,": p = ",
                                                keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   subtitle=subtitile,
                                   title = paste(method,": p = ",
                                                 keep.decmi(data.alpha$p.value),sep=""))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
  ylim.fold
  p<-p+theme(title = element_text(size=8),aspect.ratio = 1)
  if (panel.border) p<-p+theme(
    panel.border = element_rect(linetype = "solid", fill = NA,color = "black")
  )
  ggsave(paste(export_path,"/Parameter (",index,")_lineplot.pdf",sep = ""),
         width = 3,height = 3,p)
  cat("Parametric properities lineplot has been exported. Please check it.","\n")

  return(p)

}

"plotMicrochatParamMuti2Barplot"<-function(microchatParamobj,
                                         ylim.fold=1.1,
                                         xlabname=NULL,
                                         yaxis.italic=TRUE,
                                         strictmod=TRUE,
                                         panel.border=c("none","normal","all","axis"),
                                         method="anova",
                                         bar.border.size=1,
                                         seg=FALSE,
                                         comparison=my_comparisons,
                                         color_group=colorCustom(5,pal = "ywbu"),
                                         color_backgroud="grey90",
                                         axis.x.angle.adjust=FALSE,
                                         mytheme=NULL,
                                         ncol=3,
                                         spline.params=NULL,
                                         export_path="ss21/microbial parameteric analysis/liver_gene") {
  if (dev.cur() != 1) {
    dev.off()
  }
  select.index<-microchatParamobj$all.index
  paramfile.select<-microchatParamobj$paramfile.select
  export_path<-paste0(export_path,"/",paramfile.select)
  pp<-list()
  height.sel<-length(ncol_layout(length(select.index)))*3
  width.sel<-ncol_layout(length(select.index))[[1]]%>%length()*3
  for (t in select.index) {
    microchatParamStatobj<-calcMicrochatParamStat(microchatParamobj,
                                                  select.index=t,
                                                  strictmod=strictmod,
                                                  method=method,
                                                  comparison=comparison,
                                                  export_path=export_path)

    p<-plotMicrochatParam2Barplot(microchatParamStatobj,
                                  ylim.fold=ylim.fold,
                                  xlabname=xlabname,
                                  yaxis.italic=yaxis.italic,
                                  errorbar.line.add=TRUE,
                                  errorbar.point.size=0.1,
                                  y.point.adj=0.1,
                                  panel.border=panel.border,
                                  bar.border.size=bar.border.size,
                                  seg=seg,
                                  add.spline=FALSE,
                                  color_group=color_group,
                                  color_backgroud=color_backgroud,
                                  axis.x.angle.adjust=axis.x.angle.adjust,
                                  mytheme=mytheme,
                                  spline.params=spline.params,
                                  export_path=export_path)
    pp[[t]]<-p
  }

  library(patchwork)
  pm<-wrap_plots(pp,ncol=ncol)

  rownum<-length(select.index)/ncol
  if (rownum>as.integer(rownum)) nrows<-as.integer(rownum)+1 else nrows<-as.integer(rownum)

  ggsave(paste(export_path,"/Muti-parameter","_hollow_barplot.pdf",sep = ""),
         units = "cm",
         width = 7*ncol,height = 7*nrows,
         pm)
  message("Muti-parametric properities barplot has been exported. Please check it.")

  return(pm)
}

"plotMicrochatParam2Barplot"<-function(microchatParamStatobj,
                                     ylim.fold=1.1,
                                     xlabname=NULL,
                                     yaxis.italic=TRUE,
                                     errorbar.pos.adj=TRUE,
                                     errorbar.line.add=FALSE,
                                     panel.border=c("none","normal","all","axis"),
                                     seg=TRUE,
                                     add.spline=TRUE,
                                     errorbar.point.size=0.1,
                                     bar.border.size=1,
                                     y.point.adj=NULL,
                                     color_group=colorCustom(5,pal = "gygn"),
                                     color_backgroud="grey90",
                                     axis.x.angle.adjust=FALSE,
                                     mytheme=NULL,
                                     spline.params=NULL,
                                     export_path="microbial diversity analysis") {
  panel.border<-match.arg(panel.border)
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

  if (method %in% c("anova","kruskal.test")) {
    opc.data <- microchatParamStatobj$opc.data
    o.data <- opc.data[which(opc.data$index == index), ]
    gnum <- ncol(o.data)
    gnum/2
    liner <- paste(colnames(o.data[gnum/2 + 2]), " effect: ",
                   o.data[gnum/2 + 2], sep = "")
    quad <- paste(colnames(o.data[gnum/2 +3]), " effect: ",
                  o.data[gnum/2 +  3], sep = "")
    subtitile <- paste(liner, "\n",quad, sep = "")
  }

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
      geom_bar(data = data_err,aes(x = group,y = mean,color = group),
               size=bar.border.size,
               stat = "identity",position = "dodge",fill=NA,
               width = 0.7)

    if (!seg) data.seg <- calcSegmentOrder(data.alpha, y.ratio,
                                           data_poi)
    if (seg) data.seg <-calcorderx(data.alpha, y.ratio,
                                   data_poi)
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

      data.seg$y1[tt]<-y1.use
      data.seg$y2[tt]<-y2.use
      data.seg$y3[tt]<-y3.use
      data.seg$y4[tt]<-y4.use

    }

    for (tt in 1:nrow(data.seg)) {
      x1<-data.seg$order1[tt]
      x2<-data.seg$order2[tt]

      y1.use<-data.seg$y1[tt]
      y2.use<-data.seg$y2[tt]
      y3.use<-data.seg$y3[tt]
      y4.use<-data.seg$y4[tt]

      label.use=data.seg$sig[tt]

      p<-p+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x1, y = y1.use, ymin = y1.use, ymax = y1.use)+
        annotate("pointrange", color="grey50",size=errorbar.point.size,
                 x = x2, y = y2.use, ymin = y2.use, ymax = y2.use)+

        annotate("segment",color="grey50", x = x1, y = y1.use, xend = x1, yend = y3.use)+
        annotate("segment",color="grey50", x = x2, y = y2.use, xend = x2, yend = y3.use)+
        annotate("segment",color="grey50", x = x1, y = y3.use, xend = x2, yend = y3.use)+
        annotate("text",family="serif",x = (x1+x2)/2, y = y3.use+y.ratio*2/3,label=label.use)

      if (errorbar.line.add) p<-p+annotate("segment", color="grey50",x = x1-0.05, y = y1.use, xend = x1+0.05, yend = y1.use)+
        annotate("segment",color="grey50", x = x2-0.05, y = y2.use, xend = x2+0.05, yend = y2.use)

    }

    p<-p+scale_discrete_manual(values=colors,
                               aesthetics = "colour")+
      scale_discrete_manual(values=colors,
                            aesthetics = "fill")


    if (yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
            title = element_text(family = "serif",face="bold.italic", size=12))

    if (!yaxis.italic) p<-p+labs(y=index,title = method)+
      theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

      p<-ggplot() +
        geom_errorbar(data = data_err,
                      aes(x = group, y = mean, group = index,
                          ymin = mean-se, ymax = mean+se),
                      colour = "grey50",width=.25)+
        geom_bar(data = data_err,aes(x = group,y = mean,color = group),
                 stat = "identity",position = "dodge",fill=NA,
                 size=bar.border.size,
                 width = 0.7) +
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 3,color = "black",family = "serif", fontface="bold")


      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")


      if (yaxis.italic) p<-p+labs(y=index,
                                  title = paste(
                                    paste(method,": p = ",
                                          keep.decmi(data.alpha$p.value),sep=""),"\n",
                                    subtitile,sep = ""
                                  ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   title = paste(
                                     paste(method,": p = ",
                                           keep.decmi(data.alpha$p.value),sep=""),"\n",
                                     subtitile,sep = ""
                                   ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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
      if (microchatParamStatobj$param.stats$p.value>0.05) sig_label_new$alpha<-''
      data_errx<-subset(data_err, select=c(group,se))
      if ("se" %in% colnames(sig_label_new)) sig_label_new<-sig_label_new else sig_label_new<-merge(sig_label_new,data_errx,by="group")
      if(length(unique(sig_label_new$alpha))==1) sig_label_new$alpha<-""

      p<-ggplot() +geom_errorbar(data = data_err,
                                 aes(x = group, y = mean, group = index,
                                     ymin = mean-se, ymax = mean+se),
                                 colour = "grey50",width=.25)+
        geom_bar(data = data_err,aes(x = group,y = mean,color = group),
                 stat = "identity",position = "dodge",fill=NA,
                 size=bar.border.size,
                 width = 0.7) +
        geom_text(data = sig_label_new,vjust=-0.5,
                  aes(x = group,y = mean+se,label = alpha),
                  size = 3,color = "black",family = "serif", fontface="bold")

      p<-p+scale_discrete_manual(values=colors,
                                 aesthetics = "colour")+
        scale_discrete_manual(values=colors,
                              aesthetics = "fill")

      if (yaxis.italic) p<-p+labs(y=index,
                                  title = paste(
                                    paste(method,": p = ",
                                          keep.decmi(data.alpha$p.value),sep=""),"\n",
                                    subtitile,sep = ""
                                  ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold.italic",family = "serif",vjust = 1.5),
              title = element_text(face = "bold.italic",family = "serif", size=12))

      if (!yaxis.italic) p<-p+labs(y=index,
                                   title = paste(
                                     paste(method,": p = ",
                                           keep.decmi(data.alpha$p.value),sep=""),"\n",
                                     subtitile,sep = ""
                                   ))+
        theme(axis.title.y = element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
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

  ##add fitted curve for  quadratic model
  if (method %in% c("anova","kruskal.test") & !is.null(spline.params)){
    opc.thres=spline.params$opc.thres
    spline.size = spline.params$spline.size
    spline.color=spline.params$spline.color
    spline.shape =spline.params$spline.shape
    spline.type=spline.params$spline.type
    opt.interval=spline.params$opt.interval
    optimal.pt.size=spline.params$optimal.pt.size
    optimal.pt.color=spline.params$optimal.pt.color
    optimal.pt.type=spline.params$optimal.pt.type

    vline.color=spline.params$vline.color
    vline.size=spline.params$vline.size
    vline.type=spline.params$vline.type

    csal<-strsplit(o.data$Linear,"( )")[[1]][2]
    lin.num<-substr(csal,start = 2,stop =nchar(csal)-1 )%>%as.numeric()

    csa<-strsplit(o.data$Quadratic,"( )")[[1]][2]
    qua.num<-substr(csa,start = 2,stop =nchar(csa)-1 )%>%as.numeric()
    if (!is.na(lin.num)) {
    if ("Cubic" %in% colnames(o.data)) {
      csac<-strsplit(o.data$Cubic,"( )")[[1]][2]
      cub.num<-substr(csac,start = 2,stop =nchar(csac)-1 )%>%as.numeric()

      all.opc.r<-c(lin.num,qua.num,cub.num)

    } else {
      all.opc.r<-c(lin.num,qua.num)
    }

    #拟合三次函数
      if ("Cubic" %in% colnames(o.data))  if (cub.num>opc.thres & max(all.opc.r)==cub.num) {
      data_err_test<-cub_equa_fitTest(microchatParamStatobj$data_err)

      asd<-data_err_test$asd
      opt.dose<-data_err_test$opt.dose
      a<-data_err_test$a
      b<-data_err_test$b
      c<-data_err_test$c
      d<-data_err_test$d
      opt.inclusion<-opt.interval[as.integer(opt.dose)]+(opt.dose-as.integer(opt.dose))*(opt.interval[as.integer(opt.dose)+1]-opt.interval[as.integer(opt.dose)])

      l <- list(d = as.numeric(format(d, digits = 2)),
                a = as.numeric(format(a, digits = 2)),
                b = as.numeric(format(b, digits = 2)),
                c = as.numeric(format(c, digits = 2)),
                dose = as.numeric(format(opt.inclusion, digits = 2)),
                r2 = as.numeric(format(data_err_test$R2, digits = 2)))
      #eq<-substitute(italic(y) ==  d %.% (italic(x)-dose)^3 + a %.% (italic(x)-dose)^2 + b %.% (italic(x)-dose) + c~","~italic(R)^2~"="~r2,l)
      eq<-""

      p<-p + geom_line(data = data.frame(x = data_err_test$x_pred,
                                         y = data_err_test$y_pred),
                       aes(x, y),
                       color = spline.color,
                       linewidth = spline.size/2,
                       linetype=spline.type,
                       lineend = "round")+
        annotate(geom="segment",x = opt.dose, xend = opt.dose, size=vline.size,linetype=vline.type,
                 lineend = "round",
                 y =-Inf, yend = d*opt.dose^3+ a*opt.dose^2+b*opt.dose+c,color=vline.color)+
        geom_point(aes(x=opt.dose,y=d*opt.dose^3+a*opt.dose^2+b*opt.dose+c),
                   shape=optimal.pt.type,
                   size=optimal.pt.size,color=optimal.pt.color)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(d*opt.dose^3+a*opt.dose^2+b*opt.dose+c)),
                                 aes(x=opt.dose,y=y,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                 arrow = arrow(length = unit(0.01, "npc")),family="serif",size=2.5,
                                 box.padding = 1)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)*0.8),
                                 family="serif",size=2,
                                 aes(x = opt.dose, y = -Inf, label = as.character(as.expression(eq))),parse = TRUE)
    }
    #拟合二次函数
    if (qua.num>opc.thres & max(all.opc.r)==qua.num) {
      data_err_test<-qua_equa_fitTest(microchatParamStatobj$data_err)
      asd<-data_err_test$asd
      opt.dose<-data_err_test$opt.dose
      a<-data_err_test$a
      b<-data_err_test$b
      c<-data_err_test$c
      opt.inclusion<-opt.interval[as.integer(opt.dose)]+(opt.dose-as.integer(opt.dose))*(opt.interval[as.integer(opt.dose)+1]-opt.interval[as.integer(opt.dose)])

      l <- list(a = as.numeric(format(a, digits = 3)),
                b = as.numeric(format(b, digits = 3)),
                c=as.numeric(format(c, digits = 3)),
                dose = as.numeric(format(opt.inclusion, digits = 2)),
                r2 = as.numeric(format(data_err_test$R2, digits = 3)))
      #eq<-substitute(italic(y) == a %.% (italic(x)-dose)^2 + b %.% (italic(x)-dose) + c~","~italic(R)^2~"="~r2,l)
      eq<-""

      p<-p + ggalt::geom_xspline(data = asd,
                                 linetype=spline.type,lineend = "round",
                                 size = spline.size, color = spline.color,
                                 aes(x = gg, y = y_new), spline_shape = spline.shape) +
        annotate(geom="segment",x = opt.dose, xend = opt.dose, size=vline.size,linetype=vline.type,
                 lineend = "round",
                 y =-Inf, yend =  a*opt.dose^2+b*opt.dose+c,color=vline.color)+
        geom_point(aes(x=opt.dose,y=a*opt.dose^2+b*opt.dose+c),
                   shape=optimal.pt.type,
                   size=optimal.pt.size,color=optimal.pt.color)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)),
                                 aes(x=opt.dose,y=y,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                 arrow = arrow(length = unit(0.01, "npc")),family="serif",size=2.5,
                                 box.padding = 1)+
        ggrepel::geom_text_repel(data=data.frame(opt.dose=opt.dose,y=(a*opt.dose^2+b*opt.dose+c)*0.8),
                                 family="serif",size=2,
                                 aes(x = opt.dose, y = -Inf, label = as.character(as.expression(eq))),parse = TRUE)
    }
    #拟合线性模型
    if (lin.num>opc.thres & max(all.opc.r)==lin.num) {
      data_err_test<-linear_equa_fitTest(microchatParamStatobj$data_err)

      if (data_err_test$a>=0) {
        message("The fitted curve was drawn as user required. However, the slope of the fitted curve is higher than 0, so the fitted curve will not display.")
      } else {
        message("The fitted curve was drawn as user required. The slope of the fitted curve is lower than 0, so the fitted curve will be divided into two segments.")
        message("Namely, fitted curve for broken-line model is showing  !!!")

        slope<-c()
        for (t in 1:(nrow(data_err_test$asd)-1)) {
          slope<-c(slope,obtainSlope(data<-data_err_test$asd[t:(t+1),]))
        }
        slope_a<-slope[which(slope<0)]
        slope_a<-slope_a[which(slope_a==min(slope_a))]
        spot_start<-which(slope==slope_a)
        if (which(data_err_test$asd$mean==max(data_err_test$asd$mean))!=as.numeric(spot_start)) {
          spot_start<-which(data_err_test$asd$mean==max(data_err_test$asd$mean))
        } else {
          spot_start<-spot_start
        }
        if (spot_start ==2 | spot_start==3) {
          spot_end<-spot_start+1

          data_err_test1<-linear_equa_fitTest(microchatParamStatobj$data_err[1:spot_start,])
          data_err_test2<-linear_equa_fitTest(microchatParamStatobj$data_err[spot_end:nrow(microchatParamStatobj$data_err),],spot_end)

          crossdot<-getCrossCoord(data_err_test1,data_err_test2)

          if (crossdot$x>1 & crossdot$x <nrow(data_err_test$asd)) {
            asd1<-data_err_test1$asd
            a1<-data_err_test1$a
            b1<-data_err_test1$b
            asd2<-data_err_test2$asd
            asd2$gg<-(max(asd1$gg)+1):nrow(data_err_test$asd)
            a2<-data_err_test2$a
            b2<-data_err_test2$b

            asd1_a<-asd1[1:2,]
            asd1_a$gg[1]<-min(asd1$gg)-0.3
            asd1_a$y_new[1]<-a1*asd1_a$gg[1]+b1
            asd1_a$gg[2]<-crossdot$x+0.3
            asd1_a$y_new[2]<-a1*asd1_a$gg[2]+b1
            asd1<-rbind(asd1,asd1_a)

            asd2_a<-asd2[1:2,]
            asd2_a$gg[1]<-max(asd2$gg)+0.3
            asd2_a$y_new[1]<-a2*asd2_a$gg[1]+b2
            asd2_a$gg[2]<-crossdot$x-0.3
            asd2_a$y_new[2]<-a2*asd2_a$gg[2]+b2
            asd2<-rbind(asd2,asd2_a)

            opt.inclusion<-opt.interval[as.integer(crossdot$x)]+(crossdot$x-as.integer(crossdot$x))*(opt.interval[as.integer(crossdot$x)+1]-opt.interval[as.integer(crossdot$x)])

            la <- list(a = as.numeric(format(a1, digits = 2)),
                       b = as.numeric(format(b1, digits = 2)),
                       dose = as.numeric(format(opt.inclusion, digits = 2)),
                       r2 = as.numeric(format(data_err_test1$R2, digits = 2)))
            #eq1<-substitute(italic(y) ==  a %.% (italic(x)-dose) + b~","~italic(R)^2~"="~r2,la)
            lb <- list(a = as.numeric(format(a2, digits = 2)),
                       b = as.numeric(format(b2, digits = 2)),
                       dose = as.numeric(format(opt.inclusion, digits = 2)),
                       r2 = as.numeric(format(data_err_test2$R2, digits = 2)))
            #eq2<-substitute(italic(y) ==  a %.% (italic(x)-dose) + b~","~italic(R)^2~"="~r2,lb)
            eq1<-""
            eq2<-""
            p<-p + geom_line(data = asd1,
                             aes(gg, y_new),
                             color = spline.color,
                             linewidth = spline.size/2,
                             linetype=spline.type,
                             lineend = "round")+
              geom_line(data = asd2,
                        aes(gg, y_new),
                        color = spline.color,
                        linewidth = spline.size/2,
                        linetype=spline.type,
                        lineend = "round")+
              annotate(geom="segment",x = crossdot$x, xend = crossdot$x,
                       size=vline.size,linetype=vline.type,
                       lineend = "round",
                       y =-Inf, yend = a1*crossdot$x+b1,
                       color=vline.color)+
              geom_point(aes(x=crossdot$x,y=a1*crossdot$x+b1),
                         shape=optimal.pt.type,
                         size=optimal.pt.size,color=optimal.pt.color)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=crossdot$x,y=a1*crossdot$x+b1),
                                       aes(x=crossdot$x,y=-Inf,label=paste0("Optimal level: ",round(opt.inclusion,2))),
                                       arrow = arrow(length = unit(0.04, "npc")),family="serif",size=2.5,
                                       box.padding = 1)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=(crossdot$x+min(asd1$gg))/2,
                                                       y=(a1*crossdot$x+b1+0)*0.6),
                                       family="serif",size=2,arrow = arrow(length = unit(0.03, "npc")),
                                       aes(x = opt.dose, y = y, label = as.character(as.expression(eq1))),parse = TRUE)+
              ggrepel::geom_text_repel(data=data.frame(opt.dose=(crossdot$x+max(asd2$gg))/2,
                                                       y=(a2*crossdot$x+b2+0)*0.4),
                                       family="serif",size=2,arrow = arrow(length = unit(0.03, "npc")),
                                       aes(x = opt.dose, y = y, label = as.character(as.expression(eq2))),parse = TRUE)
          }
        }
      }
    }
    }
  }

  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-microchatParamStatobj$data_err$group%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p<-p+ scale_x_discrete(labels = xlabname)
  }

  p<-p+theme(aspect.ratio = 1,
             title =  element_text(size = 7,face="bold.italic"),
             axis.title = element_text(size = 6,face="bold"),
             axis.text = element_text(size = 4,face="bold"))

  if (method %in% c("anova","kruskal.test")) {
    y.offset<-1/20*max(data_poi$value)
  }

  if (method %in% c("t.test","wilcox.test")) {
    y.offset<-max(data.seg$y3)-max(data_poi$value)+y.ratio*2/3
  }

  if (panel.border=="all") p<-p+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = -Inf+1, xend = -Inf+1)+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = Inf, yend = Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = Inf+1, xend = Inf+1)

  if (panel.border=="axis") p<-p+
    annotate(geom = "segment",size=2,lineend = "round",
             x = 1-0.7/2, xend = length(unique(data_poi$group))+0.7/2,
             y = -Inf, yend = -Inf)+
    annotate(geom = "segment",size=2,lineend = "round",
             y = 0, yend = max(data_poi$value)+y.offset,
             x = -Inf+1, xend = -Inf+1)

  if (panel.border=="normal") p<-p+theme(
    panel.border = element_rect(linetype = "solid", fill = NA,color = "black",linewidth=1.5)
  )
  if (panel.border=="none") p<-p

  if (axis.x.angle.adjust) p<-p+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
  if (!is.null(mytheme)) p<-p+mytheme
  ggsave(paste(export_path,"/Parameter (",index,")_hollow_barplot.pdf",sep = ""),
         units = "cm",
         width = 21/3,
         height = 21*p$theme$aspect.ratio/3,
         p)
  cat("Parametric properities barplot has been exported. Please check it.","\n")
  cat("----------------------------------------------------------------------","\n")
  return(p)

}

"plotMicrochatParamMutiBarplot" <- function(microchatParamobj,
                                            ylim.fold=1.1,
                                            xlabname=NULL,
                                            yaxis.italic=TRUE,
                                            strictmod=TRUE,
                                            panel.border=c("none","normal","all","axis"),
                                            method="anova",
                                            seg=FALSE,
                                            comparison=my_comparisons,
                                            color_group=colorCustom(5,pal = "ywbu"),
                                            color_backgroud="grey90",
                                            axis.x.angle.adjust=FALSE,
                                            mytheme=NULL,
                                            ncol=3,
                                            spline.params=NULL,
                                            export_path="ss21/microbial parameteric analysis/liver_gene") {
  if (dev.cur() != 1) {
    dev.off()
  }
  select.index<-microchatParamobj$all.index
  paramfile.select<-microchatParamobj$paramfile.select
  export_path<-paste0(export_path,"/",paramfile.select)
  pp<-list()
  height.sel<-length(ncol_layout(length(select.index)))*3
  width.sel<-ncol_layout(length(select.index))[[1]]%>%length()*3
  for (t in select.index) {
    microchatParamStatobj<-calcMicrochatParamStat(microchatParamobj,
                                                  select.index=t,
                                                  strictmod=strictmod,
                                                  method=method,
                                                  comparison=comparison,
                                                  export_path=export_path)

    p<-plotMicrochatParamBarplot(microchatParamStatobj,
                                 ylim.fold=ylim.fold,
                                 xlabname=xlabname,
                                 yaxis.italic=yaxis.italic,
                                 errorbar.line.add=TRUE,
                                 errorbar.point.size=0.1,
                                 y.point.adj=0.1,
                                 seg=seg,
                                 add.spline=FALSE,
                                 panel.border=panel.border,
                                 color_group=color_group,
                                 color_backgroud=color_backgroud,
                                 axis.x.angle.adjust=axis.x.angle.adjust,
                                 mytheme=mytheme,
                                 spline.params=spline.params,
                                 export_path=export_path)
    pp[[t]]<-p
  }

  library(patchwork)
  pm<-wrap_plots(pp,ncol=ncol)
  rownum<-length(select.index)/ncol
  if (rownum>as.integer(rownum)) nrows<-as.integer(rownum)+1 else nrows<-as.integer(rownum)

  ggsave(paste(export_path,"/Muti-parameter","_barplot.pdf",sep = ""),
         units = "cm",
         width = 7*ncol,height = 7*nrows,
         pm)
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
    match(x,colnames(params_table)) #grep(x,colnames(params_table))
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


"calcMicrochat2Param" <- function(params) {
  if (class(params)[1]=="microchat") {
    params<-params[,1:(ncol(params)-4)]%>%data.frame()
    rownames(params)<-params$index
    params<-params[,-which(colnames(params)=="index")]
  } else {
      params<-params
      }

  params.n<-data.frame()
  times<-0
  for (i in 1:ncol(params)) {
    for (j in 1:nrow(params)) {
      params.newd<-strsplit(remove.letter(params[j,i]), "±",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
      params.newlet<-keep.letter(params[j,i])
      params.newd$letter<-params.newlet
      colnames(params.newd)<-c("mean","se","letter")
      params.newd$group<-colnames(params)[i]
      params.newd$variable<-rownames(params)[j]
      times<-times+1
      rownames(params.newd)<-times
      params.newd<-subset(params.newd,select=c(4,5,1,2,3))
      {
        params.n<-rbind(params.n,params.newd)
      }
    }
  }

  params.n$group<-factor(params.n$group,levels = unique(colnames(params)))
  params.n$variable<-factor(params.n$variable,levels = unique(rownames(params)))
  params.n<-params.n[order(params.n$variable,params.n$group),]

  params.n$mean<-params.n$mean%>%as.numeric()
  params.n$se<-params.n$se%>%as.numeric()
  return(params.n)
}
"plotMicrochatParamBoxplot1" <- function(params.n,color_bar,label.x.angle=45) {
  if (length(color_bar)!=length(unique(params.n$group))) stop("Please provide colours number equal to the groups number!!!")
  names(color_bar)<-unique(params.n$group)
  p<-ggplot(params.n,
            aes(variable,mean,fill = group))+
    geom_errorbar(data = params.n,position = position_dodge(0.9),
                  aes(x = variable, y = mean,
                      group = group, ymin = mean-se, ymax = mean+se),
                  colour = "grey50",width=.25)+
    geom_bar(aes(fill=group),color="white",
             stat = "identity",
             alpha=1,
             position = "dodge",
             width = 0.9)+
    geom_text(aes(x=variable,y=params.n$mean+params.n$se,label=letter),
              stat = "identity",
              vjust=-0.5,
              position = position_dodge(0.9),
              color="black",family="serif")+
    ylim(NA,max(params.n$mean)*1.4)+
    scale_fill_manual(values = color_bar,name="Group")+
    labs(y="Parameters")+
    theme(legend.position = "right",
          aspect.ratio = 1,
          axis.line.y.left  = element_line(arrow = arrow(length = unit(0.1, "inches"),
                                                         ends = "last", type = "open")),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = label.x.angle,vjust = 1,hjust = 1),
          text = element_text(family = "serif",face = "bold"))
  return(p)
}
"plotMicrochatParamBoxplot2" <- function(params.n,color_bar,label.x.angle=45,mytheme=NULL) {
  if (length(color_bar)!=length(unique(params.n$group))) stop("Please provide colours number equal to the groups number!!!")
  pp<-list()
  for (index in unique(params.n$variable)) {
    params.nn<-params.n[which(params.n$variable==index),]
    names(color_bar)<-unique(params.nn$group)
    p<-ggplot(params.nn,
              aes(group,mean,fill = group))+
      geom_errorbar(data = params.nn,position = position_dodge(0.9),
                    aes(x = group, y = mean,
                        group = group, ymin = mean-se, ymax = mean+se),
                    colour = "grey50",width=.25)+
      geom_bar(aes(fill=group),color="white",
               stat = "identity",
               alpha=1,
               #position = "dodge",
               width = 0.9)+
      geom_text(aes(x=group,y=mean+se,label=letter),
                #stat = "identity",
                vjust=-0.5,
                #position = position_dodge(0.9),
                color="black",family="serif")+
      scale_fill_manual(values = color_bar,name="Group")+
      ylim(NA,max(params.nn$mean)*1.2)+
      labs(y=index)+
      theme(legend.position = "right",aspect.ratio = 1,
            #axis.line.y.left  = element_line(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = label.x.angle,vjust = 1,hjust = 1),
            text = element_text(family = "serif",face = "bold"))
    if (!is.null(mytheme)) p<-p+mytheme
    pp[[index]]<-p

  }
  pp<-patchwork::wrap_plots(pp,guides = "collect")
  return(pp)
}
