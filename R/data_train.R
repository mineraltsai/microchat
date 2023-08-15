"calcMean" <- function(abunx,decim) {

  otu_rare<-abunx
  colname<-colnames(otu_rare)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)


  group_num<-length(trt_id)

  groupname<-substr(colnames(otu_rare),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otu_rare))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()



  if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otu <- lapply(
      apply(
        sapply(trt_id,function(x){
          grep(x,colnames(otu_rare))
        }),2,FUN = function(x){
          otu_rare[,x]
        }),
      function(x){
        round(rowMeans(x),decim)
      })
  } else {
    split_otu <- lapply(
      lapply(
        sapply(trt_id,function(x){
          grep(x,colnames(otu_rare))
        }),FUN = function(x){
          otu_rare[,x]
        }),
      function(x){
        round(rowMeans(x),decim)
      })
  }

  split_otu2<-split_otu%>%data.frame()
  return(split_otu2)
}

"calcSd" <- function(abunx,decim) {

  otu_rare<-abunx
  colname<-colnames(otu_rare)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)


  group_num<-length(trt_id)

  groupname<-substr(colnames(otu_rare),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otu_rare))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()



  if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otu <- lapply(
      apply(
        sapply(trt_id,function(x){
          grep(x,colnames(otu_rare))
        }),2,FUN = function(x){
          otu_rare[,x]
        }),
      function(x){
        round(rowSds(as.matrix(x)),decim)
      })
  } else {
    split_otu <- lapply(
      lapply(
        sapply(trt_id,function(x){
          grep(x,colnames(otu_rare))
        }),FUN = function(x){
          otu_rare[,x]
        }),
      function(x){
        round(rowSds(as.matrix(x)),decim)
      })
  }

  split_otu2<-split_otu%>%data.frame()
  return(split_otu2)
}

"addsample" <- function(ccc1,ccc2,decim,allnum,add.sample.size,

                        maxnum=1.15,
                        minnum=0.85) {
  ggg1<-data.frame()
  ggg2<-data.frame()
  add.sample.size=add.sample.size*3
  groupsum<-colSums(ccc1)

  for (tt in 1:allnum) {
    if (tt==1) {
      for (kk in 1:add.sample.size) {

        if (kk==1) {
          ggg<-(subset(ccc1, select=tt)+
                  round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()

          while(sum(ggg)>groupsum[tt]*maxnum | sum(ggg)<groupsum[tt]*minnum) {
            ggg<-(subset(ccc1, select=tt)+
                    round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()
          }

          ggg<-round(ggg,decim)%>%data.frame()
          ggg1<-rbind(ggg1,ggg)
        } else {
          ggg<-(subset(ccc1, select=tt)+
                  round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()

          while(sum(ggg)>groupsum[tt]*maxnum | sum(ggg)<groupsum[tt]*minnum) {
            ggg<-(subset(ccc1, select=tt)+
                    round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()
          }

          ggg<-round(ggg,decim)%>%data.frame()
          {
            ggg1<-cbind(ggg1,ggg)
          }
        }
      }
      ggg2<-rbind(ggg2,ggg1)
    } else {
      for (kk in 1:add.sample.size) {

        if (kk==1) {
          ggg<-(subset(ccc1, select=tt)+
                  round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()

          while(sum(ggg)>groupsum[tt]*maxnum | sum(ggg)<groupsum[tt]*minnum) {
            ggg<-(subset(ccc1, select=tt)+
                    round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()
          }
          ggg<-round(ggg,decim)%>%data.frame()
          ggg1<-cbind(ggg1,ggg)
        } else {
          ggg<-(subset(ccc1, select=tt)+
                  round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()

          while(sum(ggg)>groupsum[tt]*maxnum | sum(ggg)<groupsum[tt]*minnum) {
            ggg<-(subset(ccc1, select=tt)+
                    round(subset(ccc2, select=tt)*runif(1,min = -1,max = 1),decim))%>%abs()
          }
          ggg<-round(ggg,decim)%>%data.frame()
          {
            ggg1<-cbind(ggg1,ggg)
          }
        }
      }
      {
        ggg2<-cbind(ggg2,ggg1)
      }
    }
  }

  return(ggg1)
}

"calcsamplenum"<-function (otu_rare) {
  colname<-colnames(otu_rare)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)
  group_num<-length(trt_id)

  groupname<-substr(colnames(otu_rare),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otu_rare))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()
  samplenum<-collen%>%as.numeric()
  return(samplenum)
}

"calcAbun" <- function(abunx) {

  otu_rare<-abunx
  colname<-colnames(otu_rare)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)


  group_num<-length(trt_id)

  groupname<-substr(colnames(otu_rare),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otu_rare))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()



  if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otu <-
      apply(
        sapply(trt_id,function(x){
          grep(x,colnames(otu_rare))
        }),2,FUN = function(x){
          otu_rare[,x]
        })
  } else {
    split_otu <-
      lapply(
        sapply(trt_id,function(x){
          grep(x,colnames(otu_rare))
        }),FUN = function(x){
          otu_rare[,x]
        })
  }

  return(split_otu)
}

"addMicrochatsample" <- function(abunx,
                                 decim=0,
                                 addTo.sample.size=12,
                                 add.sample.size=6,
                                 maxnum=1.15,
                                 minnum=0.85) {

  otu_rare<-abunx
  colname<-colnames(otu_rare)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)
  group_num<-length(trt_id)

  groupname<-substr(colnames(otu_rare),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otu_rare))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()

  #if (collen1>=addTo.sample.size) stop("Please provide more training sample number.")
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()
  samplenum<-collen%>%as.numeric()

  ccc1<-calcMean(abunx,decim)
  ccc2<-calcSd(abunx,decim)

  #if (is.null(addTo.sample.size) & is.null(add.sample.size)) stop("Please supplement addTo.sample.size and add.sample.size.")

  #if (!is.null(addTo.sample.size) & collen1 > addTo.sample.size) stop("Please supplement more addTo.sample.size.")

  if (!is.null(addTo.sample.size) & is.null(add.sample.size)) {
    ## addTo.sample.size
    fff<-addsample(ccc1,ccc2,decim,allnum,addTo.sample.size,maxnum,minnum)

    newfff<-calcAbun(fff)
    addnum<-addTo.sample.size-collen%>%as.numeric()
    orinum<-collen%>%as.numeric()

    tkf1<-newfff[[1]][,1:addnum[1]]
    colnames(tkf1)<-paste(unique(substr(colnames(tkf1),start = 1,stop = 2)),(orinum[1]+1):(orinum[1]+addnum[1]),sep="")

    for (tk in 2:length(newfff)) {
      tkf<-newfff[[tk]][,1:addnum[tk]]
      colnames(tkf)<-
        paste(unique(substr(colnames(tkf),
                            start = 1,stop = 2)),
              (orinum[tk]+1):(orinum[tk]+addnum[tk]),sep="")

      {
        tkf1<-cbind(tkf1,tkf)
        }
    }

    abunxx<-cbind(abunx,tkf1)
    abunxxx<-calcAbun(abunxx)

    tkg1<-abunxxx[[1]]

    for (ttk in 2:length(abunxxx)) {
      tkg<-abunxxx[[ttk]]

      {
        tkg1<-cbind(tkg1,tkg)
      }
    }
  }

  if (is.null(addTo.sample.size) & !is.null(add.sample.size)) {
    ## add.sample.size
    fff<-addsample(ccc1,ccc2,decim,allnum,add.sample.size,maxnum,minnum)

    newfff<-calcAbun(fff)
    addnum<-rep(add.sample.size,length(newfff))
    orinum<-collen%>%as.numeric()

    tkf1<-newfff[[1]][,1:addnum[1]]
    colnames(tkf1)<-paste(unique(substr(colnames(tkf1),start = 1,stop = 2)),(orinum[1]+1):(orinum[1]+addnum[1]),sep="")

    for (tk in 2:length(newfff)) {
      tkf<-newfff[[tk]][,1:addnum[tk]]
      colnames(tkf)<-
        paste(unique(substr(colnames(tkf),
                            start = 1,stop = 2)),
              (orinum[tk]+1):(orinum[tk]+addnum[tk]),sep="")
      {
        tkf1<-cbind(tkf1,tkf)
        }
    }

    abunxx<-cbind(abunx,tkf1)
    abunxxx<-calcAbun(abunxx)

    tkg1<-abunxxx[[1]]

    for (ttk in 2:length(abunxxx)) {
      tkg<-abunxxx[[ttk]]

      {
        tkg1<-cbind(tkg1,tkg)
      }
    }
  }

  if (!is.null(addTo.sample.size) & !is.null(add.sample.size)) {
    ## add.sample.size
    add.sample.size=addTo.sample.size
    fff<-addsample(ccc1,ccc2,decim,allnum,add.sample.size,maxnum,minnum)

    newfff<-calcAbun(fff)
    addnum<-add.sample.size-collen%>%as.numeric()
    orinum<-collen%>%as.numeric()

    if (length(which(addnum==0))>0) {
      abunxxx<-calcAbun(abunx)

      tkg1<-abunxxx[[1]]

      for (ttk in 2:length(abunxxx)) {
        tkg<-abunxxx[[ttk]]

        {
          tkg1<-cbind(tkg1,tkg)
        }
      }


    } else {
      tkf1<-newfff[[1]][,1:addnum[1]]
      colnames(tkf1)<-paste(unique(substr(colnames(tkf1),start = 1,stop = 2)),
                            (orinum[1]+1):(orinum[1]+addnum[1]),
                            sep="")

      for (tk in 2:length(newfff)) {
        tkf<-newfff[[tk]][,1:addnum[tk]]
        colnames(tkf)<-
          paste(unique(substr(colnames(tkf),
                              start = 1,stop = 2)),
                (orinum[tk]+1):(orinum[tk]+addnum[tk]),sep="")
        {
          tkf1<-cbind(tkf1,tkf)
          }
      }

      abunxx<-cbind(abunx,tkf1)
      abunxxx<-calcAbun(abunxx)

      tkg1<-abunxxx[[1]]

      for (ttk in 2:length(abunxxx)) {
        tkg<-abunxxx[[ttk]]

        {
          tkg1<-cbind(tkg1,tkg)
        }
      }

    }


  }


  return(tkg1)
}

"addMicrochatSample"<-function(abunx,
                               corrected=TRUE,
                               robust=0.05,
                               maxnum=1.15,minnum=0.85,
                               addTo.sample.size=12,
                               add.sample.size=6){
  require(matrixStats)
  if (class(abunx)!="list") {
    if (length(unique(substr(rownames(abunx),start = 1,stop = 2)))!=1) {
      as<-strsplit(as.character(abunx[1,1]), NULL)[[1]]
      decim<-length(as)-match(".",as)
      abunx<-t(abunx)%>%data.frame()
      if (corrected) abunx<-data_corr(abunx,robust=robust,maxnum=maxnum,minnum=minnum)
      abunx<-addMicrochatsample(abunx,
                                maxnum=maxnum,minnum=minnum,
                                decim=decim,
                                addTo.sample.size=addTo.sample.size,
                                add.sample.size=add.sample.size)%>%t()%>%data.frame()
    } else {
      as<-strsplit(as.character(abunx[1,1]), NULL)[[1]]
      decim<-length(as)-match(".",as)
      if (corrected) abunx<-data_corr(abunx,robust=robust,maxnum=maxnum,minnum=minnum)
      abunx<-addMicrochatsample(abunx,
                                maxnum=maxnum,minnum=minnum,
                                decim=decim,
                                addTo.sample.size=addTo.sample.size,
                                add.sample.size=add.sample.size)

    }
  }

  if (class(abunx)=="list") {
    abunxx<-list()
    for (i in names(abunx)) {
      param_indiv<-abunx[[which(names(abunx)==i)]]
      as<-strsplit(as.character(param_indiv[1,1]), NULL)[[1]]
      decim<-length(as)-match(".",as)
      if (corrected) param_indiv<-data_corr(param_indiv,robust=robust,maxnum=maxnum,minnum=minnum)
      param_indiv<-t(param_indiv)%>%data.frame()
      param_indiv<-param_indiv%>%as.matrix()
      param_indiv[is.na(param_indiv)]<-0
      param_indiv<-param_indiv%>%data.frame()
      param_indiv<-addMicrochatsample(param_indiv,
                                      maxnum=maxnum,minnum=minnum,
                                      decim=decim,
                                      addTo.sample.size=addTo.sample.size,
                                      add.sample.size=add.sample.size)%>%t()%>%data.frame()

      abunxx[[which(names(abunx)==i)]]<-param_indiv
      names(abunxx)[which(names(abunx)==i)]<-i
    }
    abunx<-abunxx
  }
  return(abunx)
}



"data_corr" <- function(param_table,robust=0.05,maxnum=1.25,minnum=0.75) {
  param_table<-tibble::rownames_to_column(param_table,var = "sample")
  param_table$group<-substr(param_table$sample,start = 1,stop = 2)
  param_table<-subset(param_table, select = c(which(colnames(param_table)=="group"),
                                              which(colnames(param_table)!="group")))


  param_table$group<-factor(param_table$group,levels = unique(param_table$group))

  ###output statistics table
  all.index<-colnames(param_table)[3:length(colnames(param_table))]

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

  single.index.paramtabxx<-data.frame()

  for (tt in 3:ncol(param_table)) {
    if (tt==3) {
      single.index.statab<-summ_fcbv.all[which(summ_fcbv.all$index == colnames(param_table)[tt]),]
      single.index.paramtab<-subset(param_table,select = c(1,2,tt))
      for (tk in 1:nrow(single.index.statab)) {
        single.index.paramtabx<-single.index.paramtab[which(single.index.paramtab$group %in% single.index.statab$group[tk]),]
        mean.select<-single.index.statab$mean[tk]
        mean.select.up<-mean.select*maxnum
        mean.select.down<-mean.select*minnum
        for (ttk in 1:nrow(single.index.paramtabx)) {
          if (single.index.paramtabx[ttk,3]>mean.select.up | single.index.paramtabx[ttk,3]<mean.select.down ) {
            single.index.paramtabx[ttk,3]<-rnorm(1,mean = mean.select,sd=mean.select*robust)
          } else {
            single.index.paramtabx[ttk,3]<-single.index.paramtabx[ttk,3]
          }

        }

        {single.index.paramtabxx<-rbind(single.index.paramtabxx,single.index.paramtabx)}

      }
    } else {
      single.index.statab<-summ_fcbv.all[which(summ_fcbv.all$index == colnames(param_table)[tt]),]
      single.index.paramtab<-subset(param_table,select = c(1,2,tt))
      single.index.paramtabxxx<-data.frame()
      for (tk in 1:nrow(single.index.statab)) {
        single.index.paramtabx<-single.index.paramtab[which(single.index.paramtab$group %in% single.index.statab$group[tk]),]
        mean.select<-single.index.statab$mean[tk]
        mean.select.up<-mean.select*maxnum
        mean.select.down<-mean.select*minnum
        for (ttk in 1:nrow(single.index.paramtabx)) {
          if (single.index.paramtabx[ttk,3]>mean.select.up | single.index.paramtabx[ttk,3]<mean.select.down ) {
            single.index.paramtabx[ttk,3]<-rnorm(1,mean = mean.select,sd=mean.select*robust)
          } else {
            single.index.paramtabx[ttk,3]<-single.index.paramtabx[ttk,3]
          }

        }
        {single.index.paramtabxxx<-rbind(single.index.paramtabxxx,single.index.paramtabx)}
        single.index.paramtabxxxx<-single.index.paramtabxxx
      }
      {single.index.paramtabxx<-cbind(single.index.paramtabxx,single.index.paramtabxxxx)}
    }
  }

  single.index.paramtabxx
  rm.col<-which(colnames(single.index.paramtabxx)=="group"|colnames(single.index.paramtabxx)=="sample")
  rm.colx<-rm.col[3:length(rm.col)]
  single.index.paramtabxxs<-single.index.paramtabxx[,-rm.colx]
  single.index.paramtabxxs<-tibble::column_to_rownames(single.index.paramtabxxs,var="sample")

  single.index.paramtabxxs<-subset(single.index.paramtabxxs,select = -group)
  return(single.index.paramtabxxs)
}
