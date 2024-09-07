
"schur" <- function(a,b) {
  return(a*b)
}

"sourcetrack_test" <- function(abun,sink_group = "ml",
                               source_group="ct",
                               type=c("times","common")) {

  type<-match.arg(type)
  sink_group = sink_group

  if (sink_group != "all") {
    if (type=="times") {
      all_group<-unique(substr(colnames(abun),start = 1,stop = 2))
      source_group<-all_group[1:(match(sink_group,all_group)-1)]

      stdata<-sourcetracker_file(abun,source_group = source_group,
                                 sink_group = sink_group)

      metadata<-stdata$metadata
      otus<-stdata$otus

      data_st<-sourcetrack(metadata,otus)

      FEAST_output<-data_st
      colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号

      head(FEAST_output)
    }

    if (type=="common") {
      all_group<-unique(substr(colnames(abun),start = 1,stop = 2))


      stdata<-sourcetracker_file(abun,source_group = source_group,
                                 sink_group = sink_group)

      metadata<-stdata$metadata
      otus<-stdata$otus

      data_st<-sourcetrack(metadata,otus)

      FEAST_output<-data_st
      colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号

      head(FEAST_output)
    }
  }

  if (sink_group == "all") {
    if (type=="times") {

      all_group<-unique(substr(colnames(abun),start = 1,stop = 2))

      feastdata<-data.frame()
      for (tar_group in all_group[2:length(all_group)]) {
        sink_group1=tar_group

        source_group<-all_group[1:(match(sink_group1,all_group)-1)]

        stdata<-sourcetracker_file(abun,source_group = source_group,
                                   sink_group = sink_group1)

        metadata<-stdata$metadata
        otus<-stdata$otus

        data_st<-sourcetrack(metadata,otus)

        FEAST_output<-data_st
        colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号
        FEAST_output<-FEAST_output%>%rownames_to_column(var = "row")

        head(FEAST_output)
        {
          feastdata<-rbind(feastdata,FEAST_output)
        }
        head(feastdata)
      }

      feastdata1<-feastdata%>%group_by(row)%>%
        summarise_all(sum)
      feastdata1<-feastdata1%>%column_to_rownames(var = "row")
      colsumfeastdata1<-colSums(feastdata1)%>%data.frame()

      feastdata2<-feastdata1
      for (tt in 1:length(colnames(feastdata1))) {
        feastdata2[,tt]<-feastdata1[,tt]/colsumfeastdata1$.[tt]
      }

      groupnum<-match(all_group,rownames(feastdata2))
      groupnum1<-match(rownames(feastdata2),all_group)
      groupnum1[is.na(groupnum1)]<-max(groupnum1[!is.na(groupnum1)])+1
      groupnum[is.na(groupnum)]<-setdiff(groupnum1,groupnum)

      feastdata3<-feastdata2[groupnum,]

      FEAST_output<-feastdata3
    }

    if (type=="common") {
      all_group<-unique(substr(colnames(abun),start = 1,stop = 2))
      samplenum<-length(colnames(abun))

      trid_id<-substr(colnames(abun),start = 1,stop = 2)%>%unique()
      ttk<-sapply(trid_id, function(y){grep(y,colnames(abun))})

      sample.max.size<-lapply(ttk,function(x){length(x)})%>%as.numeric()%>%max()

      feastdata<-data.frame()

      for (tar_group in all_group[2:length(all_group)]) {
        sink_group1=tar_group

        stdata<-sourcetracker_file(abun,source_group = source_group,
                                   sink_group = sink_group1)

        metadata<-stdata$metadata
        otus<-stdata$otus

        data_st<-sourcetrack(metadata,otus)

        FEAST_output<-data_st
        colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号
        FEAST_output<-FEAST_output%>%rownames_to_column(var = "row")

        sample.sie<-length(colnames(FEAST_output))-1
        sample.sie.error<-sample.max.size-sample.sie

       if (sample.sie.error!=0) {
         data.error<-data.frame(
          row=rep(NA,sample.sie.error),
          unknown=rep(NA,sample.sie.error)
       )%>%t()%>%data.frame()
         colnames(data.error)<-paste("repeat_",(sample.sie+1):sample.max.size,sep = "")
         FEAST_output<-cbind(FEAST_output,data.error)
       }

        head(FEAST_output)


        {
          feastdata<-rbind(feastdata,FEAST_output)
        }
        head(feastdata)
      }

      feastdata_data1<-rowMeans(feastdata[,2:length(colnames(feastdata))],na.rm = TRUE)%>%data.frame()
      feastdata_data1[,"row"]<-feastdata$row
      colnames(feastdata_data1)[1]<-"value"

      feastdata_data1[,"sink"]<-rep(all_group[2:length(all_group)],each=2)
      FEAST_output<-feastdata_data1
      head(FEAST_output)
    }
  }

  return(list(FEAST_output=FEAST_output,sink_group=sink_group,source_group=source_group))
}


"sourcetracker_file" <- function(sink_group = "M3",
                                 abun,source_group = c("M0","M3","M4")) {

  group_all<-substr(colnames(abun),start = 1,stop = 2)

  source_num<-length(source_group)

  sink_abun<-abun[,which(group_all %in% sink_group)]
  source_abun<-abun[,which(group_all %in% source_group)]

  colnum<-ncol(sink_abun)+ncol(sink_abun)*ncol(source_abun)

  metadata<-data.frame(SampleID=1:colnum,
                       Env=1:colnum,
                       SourceSink=1:colnum,
                       id=1:colnum)

  metadata$Env[1:ncol(sink_abun)]<-sink_group
  metadata$Env[(ncol(sink_abun)+1):colnum]<-rep(rep(source_group,each=ncol(sink_abun)),ncol(source_abun)/source_num)

  metadata$SourceSink[which(metadata$Env == sink_group)] <- "Sink"
  metadata$SourceSink[which(metadata$Env != sink_group)] <- "Source"
  metadata$id[(ncol(sink_abun)+1):colnum]<-rep(1:ncol(sink_abun),each=ncol(sink_abun)*source_num)

  abun_all<-cbind(sink_abun,rep(source_abun,ncol(sink_abun)))

  colnames(abun_all)<-paste("FM",1:colnum,sep = "")

  metadata$SampleID<-colnames(abun_all)
  otus<-abun_all
  metadata<-tibble::column_to_rownames(metadata,var = "SampleID")


  return(list(otus=otus,metadata=metadata))
}



"calc.st" <- function(it,Proportions_est,Ids,different_sources_flag=1,metadata,otus,envs,EM_iterations) {

  # Extract the source environments and source/sink indices
  if(different_sources_flag == 1){

    train.ix <- which(metadata$SourceSink=='Source' & metadata$id == Ids[it])
    test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])

  } else{

    train.ix <- which(metadata$SourceSink=='Source')
    test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])
  }

  num_sources <- length(train.ix)
  COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))
  #Can be adjusted by the user

  # Define sources and sinks
  sources <- as.matrix(microchat::rarefy(otus[train.ix,], COVERAGE))
  sinks <- as.matrix(microchat::rarefy(t(as.matrix(otus[test.ix,])), COVERAGE))
  print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
  print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
  print(paste("The sink is:", envs[test.ix]))

  # Estimate source proportions for each sink

  FEAST_output<-microchat::FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
  Proportions_est <- FEAST_output$data_prop[,1]


  if(length(unique(as.character(envs[train.ix])))<length(as.character(envs[train.ix]))) {
    feastdata<-list()
    for (k in 1:length(unique(as.character(envs[train.ix])))) {
      samplename<-unique(as.character(envs[train.ix]))[k]
      pos<-which(as.character(envs[train.ix])==samplename)

      feastdata1<-Proportions_est[pos]%>%sum()
      {
        feastdata<-rbind(feastdata,feastdata1)
      }
    }

    feastdata<-feastdata%>%as.numeric()

    Proportions_est<-c(feastdata,
                             Proportions_est[length(Proportions_est)])
    names(Proportions_est) <- c(unique(as.character(envs[train.ix])), "unknown")
  }

  print("Source mixing proportions")
  print(Proportions_est)

  return(Proportions_est)
}


"sourcetrack" <- function(metadata,otus,nworker=4) {
  require(parallel)
  EM_iterations = 1000 #default value
  nworker=nworker
  ##if you use different sources for each sink, different_sources_flag = 1, otherwise = 0
  different_sources_flag = 1

  otus <- t(as.matrix(otus))

  # Extract only those samples in common between the two tables
  common.sample.ids <- intersect(rownames(metadata), rownames(otus))
  otus <- otus[common.sample.ids,]
  metadata <- metadata[common.sample.ids,]
  # Double-check that the mapping file and otu table
  # had overlapping samples
  if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                     'between the metadata file and data table')
    stop(message)
  }


  if(different_sources_flag == 0){

    metadata$id[metadata$SourceSink == 'Source'] = NA
    metadata$id[metadata$SourceSink == 'Sink'] = c(1:length(which(metadata$SourceSink == 'Sink')))
  }


  envs <- metadata$Env
  Ids <- na.omit(unique(metadata$id))
  Proportions_est <- list()

  if (nworker!=1) {
    c1 <- try(parallel::makeCluster(nworker, type = "PSOCK"))
    if (class(c1)[1] == "try-error") {
      c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))
    }
    if (class(c1)[1] == "try-error") {
      c1 <- try(parallel::makeCluster(nworker, setup_strategy = "sequential"))
    }
    message("Now randomizing by parallel computing. Begin at ",
            date(), ". Please wait...")
    clusterExport(c1, "%>%")
    Proportions_est = parallel::parLapply(c1, 1:length(Ids), calc.st,
                                          Proportions_est,Ids,
                                          different_sources_flag=1,
                                          metadata=metadata,
                                          otus=otus,
                                          envs=envs,
                                          EM_iterations=1000)

    parallel::stopCluster(c1)
  }

  ss<-lapply(Proportions_est, function(x){
   data.frame(x)

  })
  rowna<-data.frame()
  for (tt in 1:length(ss)) {
    sss<-ss[[tt]]
    rown<-rownames(sss)
    rown<-rown%>%data.frame()
    {
      rowna<-rbind(rowna,rown)
    }
  }
  unionrowname<-unique(rowna$.)


  for (k in 1:length(ss)) {
    datapro<-ss[[k]]
    addrownames<-setdiff(unionrowname,rownames(datapro))
    datapro[paste(addrownames),]<-0
    #datapro<-datapro[order(rownames(datapro)),]%>%data.frame()
    rownames(datapro)<-unionrowname
    ss[[k]]<-datapro
  }

  Proportions_est<-ss
  #输出计算结果

  FEAST_output = as.data.frame(Proportions_est)
  return(FEAST_output)
}

"change_C"<-function(newcov, X){

  X=t(as.matrix(X))
  idx = 1:dim(X)[2]

  if(sum(X) > newcov){

    while(sum(X) > newcov){
      greaterone = X > 1
      samps = 20
      if(samps > length(X[greaterone]))
        samps = length(X[greaterone])
      changeidx = sample(idx[greaterone], samps, replace = F)
      X[changeidx] = X[changeidx] - 1
    }

  }

  if(sum(X) < newcov){

    while(sum(X) < newcov){
      greaterone = X > 1
      samps = 100
      if(samps > length(X[greaterone]))
        samps = length(X[greaterone])
      changeidx = sample(idx[greaterone], samps, replace = F)
      X[changeidx] = X[changeidx] + 1
    }

  }

  return(X)
}

"rarefy" <- function(x,maxdepth){


  if(is.null(maxdepth)) x<-x

  x <- matrix(x,nrow=nrow(x))
  nr <- nrow(x)
  nc <- ncol(x)

  for(i in 1:nrow(x)){
    if(sum(x[i,]) > maxdepth){
      prev.warn <- options()$warn
      options(warn=-1)
      s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
      options(warn=prev.warn)
      x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
    }
  }
  return(x)
}

"jsdmatrix" <- function(x){
  d <- matrix(0,nrow=nrow(x),ncol=nrow(x))
  for(i in 1:(nrow(x)-1)){
    for(j in (i+1):nrow(x)){
      d[i,j] <- jsd(x[i,], x[j,])
      d[j,i] <- d[i,j]
    }
  }
  return(d)
}

"jsd" <- function(p,q){
  m <- (p + q)/2
  return((kld(p,m) + kld(q,m))/2)
}

"kld" <- function(p,q){
  nonzero <- p>0 & q>0
  return(sum(p[nonzero] * log2(p[nonzero]/q[nonzero])))
}

"h"<-function(x) {y <- x[x > 0]; -sum(y * log(y))};
"mult_JSD" <- function(p,q) {h(q %*% p) - q %*% apply(p, 1, h)}

"retrands"<-function(V){
  toret<-unlist(lapply(c(V), function(x) runif(1, x+1e-12, x+1e-09)))
  return(toret)
}

"getR2"<-function(x,y){
  return((cor(x,y))^2)
}

"EE"<-function(alphas, sources){
  nums<-(sapply(1:length(alphas), function(n) Reduce("+", crossprod(as.numeric(alphas[n]),as.numeric(sources[[n]])))))
  denom<-(Reduce("+", nums))
  return(nums/denom)
}

"A"<-function(alph, XO, raos){
  tmp<-crossprod(alph, XO/raos)
  tmp<-rapply(list(tmp), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  tmp<-Reduce("+",unlist(tmp))
  return(tmp)
}

"M"<-function(alphas, sources, sink, observed){

  newalphs<-c()
  rel_sink <-sink/sum(sink)

  if(sum(sources[[1]]) > 1){

    sources <-lapply(sources, function(x) x/(sum(colSums(x))))
  }


  LOs<-lapply(sources, schur, b=rel_sink)
  BOs<-t(mapply(crossprod, x=sources, y=alphas))
  BOs<-split(BOs, seq(nrow(BOs)))
  BOs<-lapply(BOs, as.matrix)
  BOs<-lapply(BOs, t)
  num_list <- list()
  source_new <- list()


  for(i in 1:length(sources)){
    num <- c()
    denom <- c()
    num<-crossprod(alphas[i], (LOs[[i]]/(Reduce("+", BOs))))
    num<-rapply(list(num), f=function(x) ifelse(is.nan(x),0,x), how="replace" ) #replace na with zero
    num_list[[i]]<- num[[1]][1,] + observed[[i]][1,]

    denom <- Reduce("+",unlist(num_list[[i]]))
    source_new[[i]] <- num_list[[i]]/denom
    source_new[[i]][is.na(source_new[[i]])] = 0
  }

  sources = source_new

  newalphs<-c()
  #sink<-as.matrix(sink); #src1<-as.matrix(sources[[1]]); src2<-as.matrix(sources[[2]])
  sources<-lapply(sources, t)
  XOs<-lapply(sources,schur, b=rel_sink)
  AOs<-t(mapply(crossprod, x=sources, y=alphas))
  AOs<-split(AOs, seq(nrow(AOs)))
  AOs<-lapply(AOs, as.matrix)
  AOs<-lapply(AOs, t)
  newAs<-c()
  for(i in 1:length(sources)){
    newA<-crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA<-rapply(list(newA), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    newA<-Reduce("+",unlist(newA))
    newAs<-c(newAs, newA)
  }
  tot<-sum(newAs)
  Results <- list (new_alpha = newAs/(tot), new_sources = sources)
  return(Results)
}

"do_EM"<-function(alphas, sources, observed, sink, iterations){

  curalphas<-alphas
  newalphas<-alphas
  m_guesses<-c(alphas[1])
  for(itr in 1:iterations){

    curalphas<-EE(newalphas, sources)
    tmp <- M(alphas = curalphas, sources = sources, sink = sink, observed = observed)
    newalphas <- tmp$new_alpha
    sources <- tmp$new_sources

    m_guesses<-c(m_guesses, newalphas[1])
    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break

  }
  toret<-c(newalphas)
  results <- list(toret = toret, sources = sources)

  return(results)
}

"M_basic" <-function(alphas, sources, sink){
  newalphs<-c()
  XOs<-lapply(sources,schur, b=sink)
  AOs<-t(mapply(crossprod, x=sources, y=alphas))
  AOs<-split(AOs, seq(nrow(AOs)))
  AOs<-lapply(AOs, as.matrix)
  AOs<-lapply(AOs, t)
  newAs<-c()
  for(i in 1:length(sources)){
    newA<-crossprod(alphas[i], (XOs[[i]]/(Reduce("+", AOs))))
    newA<-rapply(list(newA), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
    newA<-Reduce("+",unlist(newA))
    newAs<-c(newAs, newA)
  }
  tot<-sum(newAs)
  return(newAs/(tot))
}

"do_EM_basic"<-function(alphas, sources, sink, iterations){
  curalphas<-alphas
  newalphas<-alphas
  m_guesses<-c(alphas[1])
  for(itr in 1:iterations){
    curalphas<-EE(newalphas, sources)
    newalphas<-M_basic(curalphas, sources, sink)
    m_guesses<-c(m_guesses, newalphas[1])

    if(abs(m_guesses[length(m_guesses)]-m_guesses[length(m_guesses)-1])<=10^-6) break
  }
  toret<-c(newalphas)
  return(toret)
}

"source_process_nounknown" <- function(train, envs, rarefaction_depth=1000){

  train <- as.matrix(train)

  # enforce integer data
  if(sum(as.integer(train) != as.numeric(train)) > 0){
    stop('Data must be integral. Consider using "ceiling(datatable)" or ceiling(1000*datatable) to convert floating-point data to integers.')
  }
  envs <- factor(envs)
  train.envs <- sort(unique(levels(envs)))

  # rarefy samples above maxdepth if requested
  if(!is.null(rarefaction_depth) && rarefaction_depth > 0) train <- rarefy(train, rarefaction_depth)

  # get source environment counts
  # sources is nenvs X ntaxa
  X <- t(sapply(split(data.frame(train), envs), colSums))

  rownames(X) <- c(train.envs)
  X <- t(as.matrix(X))

  return(X)
}

"read_pseudo_data" <-function(dataset){
  path_to_data<-"../data/"
  if(dataset=="DA"){
    df<-read.table(paste0(path_to_data,"DA_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else if(dataset=="DB"){
    df<-read.table(paste0(path_to_data,"DB_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else if (dataset=="F4"){
    df<-read.table(paste0(path_to_data,"F4_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])
  }else{
    df<-read.table(paste0(path_to_data,"M3_99_T_d10000_date_nan.txt"), fill = NA)
    return(df[complete.cases(df),])}
}

"create_m" <- function(num_sources, n, EPSILON){


  if( n == 1 ){

    index = sample(c(1:num_sources), 1)
    m_1 = runif(min = 0.6, max = 0.9, n = 1)
    resid = 1-m_1
    other_ms = resid/(num_sources-1)
    m = rep(NA, num_sources)
    m[index] = c(m_1)
    m[is.na(m)] = other_ms

  }


  if( n == 2 ){

    index = sample(c(1:num_sources), 2)
    m_1 = runif(min = 0.1, max = 0.2, n = 1)
    m_2 = runif(min = 0.4, max = 0.5, n = 1)
    resid = 1-(m_1+m_2)
    other_ms = resid/(num_sources-2)
    m = rep(NA, num_sources)
    m[index] = c(m_1, m_2)
    m[is.na(m)] = other_ms

  }


  if( n == 3 ){

    index = sample(c(1:num_sources), 3)
    m_1 = runif(min = 0.1, max = 0.5, n = 1)
    m_2 = runif(min = 0.2, max = 0.25, n = 1)
    m_3 = runif(min = 0.1, max = 0.15, n = 1)
    resid = 1-(m_1+m_2+m_3)
    other_ms = runif(min = 0.001, max = resid/(num_sources-3), n = (num_sources-3))
    m = rep(NA, num_sources)
    m[index] = c(m_1, m_2, m_3)
    m[is.na(m)] = other_ms
    m = m/sum(m)

  }
  subsum = 0
  idx = 1:length(m)

  while ((subsum+0.001) < EPSILON){
    tosub = EPSILON - subsum
    tosub = tosub / (num_sources+1)
    mask = m > tosub
    m[mask] = m[mask] - tosub
    subsum = subsum + length(m[mask]) * tosub

  }
  m = c(m,(EPSILON))

  # sum(m)
  return(m)

}

"unknown_initialize" <- function(sources, sink, n_sources){

  unknown_source = rep(0, length(sink))
  sum_sources = apply(sources, 2, sum)

  unknown_source = c()

  for(j in 1:length(sum_sources)){

    unknown_source[j] = max(sink[j]-sum_sources[j], 0)

  }



  return(unknown_source)

}

"unknown_initialize_1" <- function(sources, sink, n_sources){

  unknown_source = rep(0, length(sink))
  sources_sum = apply(sources, 2 ,sum)


  unknown_source = c()

  for(j in 1:length(sources_sum)){

    unknown_source[j] = max(sink[j]-sources_sum[j], 0)

  }

  #Select the cor OTUs
  ind_cor = list()
  ind_known_source_abun = c()
  ind_cor_all = which(sources[1,] > 0)

  counter = matrix(0, ncol = dim(sources)[2], nrow =  dim(sources)[1])


  for(j in 1:n_sources){

    ind_cor[[j]] = which(sources[j,] > 0)

    for(k in 1:length(sources[j,])){

      if(sources[j,k] > 0){

        counter[j,k] = counter[j,k]+1
      }


    }

  }

  OTU_present_absent = apply(counter, 2, sum)
  ind_cor_all = which(OTU_present_absent >= round(n_sources*0.8))

  if(length(ind_cor_all) > 1){

    cor_abundance = round(apply(sources[,ind_cor_all], 2, median)/2) #take the min abundnace of the 'cor'
    unknown_source[ind_cor_all] = cor_abundance

  }


  #keep the sink abundance where there is no known source
  ind_no_known_source_abun = which(sources_sum == 0)

  for(j in 1:length(ind_no_known_source_abun)){

    # unknown_source[ind_no_known_source_abun[j]] = max(runif(n = 1, min = 1, max = 100), sink[ind_no_known_source_abun[j]])
    unknown_source[ind_no_known_source_abun[j]] = max((sink[ind_no_known_source_abun[j]] - rpois(n = 1, lambda = 0.5)), 0)

  }



  return(unknown_source)

}

"unknown__initialize_1" <- function(sources, sink, n_sources){

  unknown_source = rep(0, length(sink))

  #zero all the OTUs with at least 1 known source
  sources_sum = apply(sources, 2 ,sum)
  ind_known_source_abun = which(sources_sum > 0)
  unknown_source[ind_known_source_abun] = 0


  #Select the cor OTUs
  ind_cor = list()
  ind_known_source_abun = c()
  ind_cor_all = which(sources[1,] > 0)

  counter = matrix(0, ncol = dim(sources)[2], nrow =  dim(sources)[1])


  for(j in 1:n_sources){

    ind_cor[[j]] = which(sources[j,] > 0)

    for(k in 1:length(sources[j,])){

      if(sources[j,k] > 0){

        counter[j,k] = counter[j,k]+1
      }


    }

  }

  OTU_present_absent = apply(counter, 2, sum)
  ind_cor_all = which(OTU_present_absent >= round(n_sources*0.8))

  if(length(ind_cor_all) > 1){

    cor_abundance = apply(sources[,ind_cor_all], 2, median) #take the median abundnace of the 'cor'
    unknown_source[ind_cor_all] = cor_abundance

  }



  #keep the sink abundance where there is no known source
  ind_no_known_source_abun = which(sources_sum == 0)

  for(j in 1:length(ind_no_known_source_abun)){

    unknown_source[ind_no_known_source_abun[j]] = max( round(sink[ind_no_known_source_abun[j]]+ rnorm(n = length(sink[ind_no_known_source_abun[j]]))), 0)

  }



  return(unknown_source)

}


"FEAST" <- function(source = sources_data, sinks = sinks, em_itr = 1000, env = rownames(sources_data), include_epsilon = T, COVERAGE,
                    unknown_initialize_flag = 0){


  tmp = source
  test_zeros = apply(tmp, 1, sum)
  ind_to_use = as.numeric(which(test_zeros > 0))
  ind_zero = as.numeric(which(test_zeros == 0))

  source = tmp[ind_to_use,]
  sinks = sinks



  #####adding support for multiple sources#####
  totalsource<-source
  totalsource<-as.matrix(totalsource)
  sources <- split(totalsource, seq(nrow(totalsource)))
  sources<-lapply(sources, as.matrix)
  dists<-lapply(sources, function(x) x/(sum(colSums(x))))
  totaldist<-t(Reduce("cbind", dists))
  sinks<-matrix(sinks, nrow = 1, ncol = dim(totalsource)[2])

  num_sources = dim(source)[1]
  envs_simulation = c(1:(num_sources))

  source_old = source
  totalsource_old = totalsource

  source_old=lapply(source_old,t)
  source_old<- split(totalsource_old, seq(nrow(totalsource_old)))
  source_old<-lapply(source_old, as.matrix)

  #Creating the unknown source per mixing iteration
  if(include_epsilon == TRUE){

    ##Adding the initial value of the unknown source for CLS and EM
    source_2 = list()
    totalsource_2 = matrix(NA, ncol = dim(totalsource_old)[2], nrow = ( dim(totalsource_old)[1] + 1))

    for(j in 1:num_sources){

      source_2[[j]] = source_old[[j]]
      totalsource_2[j,] = totalsource_old[j,]
    }

    #create unknown for each sink i

    sinks_rarefy = rarefy(matrix(sinks, nrow = 1), maxdepth = apply(totalsource_old, 1, sum)[1]) #make

    if(unknown_initialize_flag == 1)
      unknown_source = unknown_initialize_1(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                            n_sources = num_sources)


    if(unknown_initialize_flag == 0)
      unknown_source = unknown_initialize(sources = totalsource[c(1:num_sources),], sink = as.numeric(sinks),
                                          n_sources = num_sources)

    # unknown_source = unknown_source_1 + rpois(n = length(sinks), lambda = 0.5)

    unknown_source_rarefy = rarefy(matrix(unknown_source, nrow = 1), maxdepth = COVERAGE)
    source_2[[j+1]] = t(unknown_source_rarefy)
    totalsource_2[(j+1),] = t(unknown_source_rarefy)
    totalsource = totalsource_2

    source=lapply(source_2,t)
    # totalsource <- rarefy(x = totalsource, maxdepth = COVERAGE)
    source<- split(totalsource, seq(nrow(totalsource_2)))
    source<-lapply(source_2, as.matrix)

    envs_simulation <- c(1:(num_sources+1))

  }


  samps <- source
  samps<-lapply(samps, t)

  observed_samps <- samps
  observed_samps[[(num_sources + 1)]] = t(rep(0, dim(samps[[1]])[2]))


  #Calculate JSD value
  # x <- totalsource[c(1:num_sources),]
  # JSDMatrix <- jsdmatrix(x)
  # JSDMatrix <- JSDMatrix/COVERAGE
  # JS = mean(JSDMatrix[-which(JSDMatrix == 0)])
  # js_values = append(js_values, JS)
  # print(js_values)

  initalphs<-runif(num_sources+1, 0.0, 1.0)
  initalphs=initalphs/Reduce("+", initalphs)
  sink_em = as.matrix(sinks)
  pred_em<-do_EM_basic(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr)

  tmp<-do_EM(alphas=initalphs, sources=samps, sink=sink_em, iterations=em_itr, observed=observed_samps)
  pred_emnoise = tmp$toret

  k = 1
  pred_emnoise_all = c()
  pred_em_all = c()

  for(j in 1:length(env)){

    if(j %in% ind_to_use){

      pred_emnoise_all[j] = pred_emnoise[k]
      pred_em_all[j] = pred_em[k]
      k = k+1

    }

    else{

      pred_emnoise_all[j] = 0
      pred_em_all[j] = 0
    }

  }

  pred_emnoise_all[j+1] = pred_emnoise[k]
  pred_em_all[j+1] = pred_em[k]



  names(pred_emnoise_all) = c(env,"unknown")
  names(pred_em_all) = c(env,"unknown")


  Results = list(unknown_source = unknown_source, unknown_source_rarefy = unknown_source_rarefy,
                 data_prop = data.frame(pred_emnoise_all,pred_em_all))
  return(Results)

}




"calcMicrochatSourceTimes" <- function(submchat,
                                       nworker=4,
                                       sink_group = "ml",
                                       source_group="ct") {


  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  sink_group = sink_group
  abun<-submchat$otu_table

  if (sink_group != "all") {

    all_group<-unique(substr(colnames(abun),start = 1,stop = 2))
    source_group<-all_group[1:(match(sink_group,all_group)-1)]

    stdata<-sourcetracker_file(abun,source_group = source_group,
                               sink_group = sink_group)

    metadata<-stdata$metadata
    otus<-stdata$otus

    data_st<-sourcetrack(metadata,otus,nworker=nworker)

    FEAST_output<-data_st
    colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号

    head(FEAST_output)

  }

  if (sink_group == "all") {

    all_group<-unique(substr(colnames(abun),start = 1,stop = 2))

    samplenum<-length(colnames(abun))

    trid_id<-substr(colnames(abun),start = 1,stop = 2)%>%unique()
    ttk<-sapply(trid_id, function(y){grep(y,colnames(abun))})

    sample.max.size<-lapply(ttk,function(x){length(x)})%>%as.numeric()%>%max()


    feastdata<-data.frame()
    for (tar_group in all_group[2:length(all_group)]) {
      sink_group1=tar_group

      source_group<-all_group[1:(match(sink_group1,all_group)-1)]

      stdata<-sourcetracker_file(abun,source_group = source_group,
                                 sink_group = sink_group1)

      metadata<-stdata$metadata
      otus<-stdata$otus

      data_st<-sourcetrack(metadata,otus,nworker=4)

      FEAST_output<-data_st
      colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号
      FEAST_output<-FEAST_output%>%rownames_to_column(var = "row")

      sample.sie<-length(colnames(FEAST_output))-1
      sample.sie.error<-sample.max.size-sample.sie

      if (sample.sie.error!=0) {
        rowmean.use<-rowMeans(FEAST_output[,2:ncol(FEAST_output)])
        row.na<-length(rowmean.use)
        data.error<-rowmean.use%>%data.frame()
        rownames(data.error)<-rownames(FEAST_output)
        colnames(data.error)<-paste("repeat_",(sample.sie+1):sample.max.size,sep = "")
        FEAST_output<-cbind(FEAST_output,data.error)
      }


      head(FEAST_output)
      {
        feastdata<-rbind(feastdata,FEAST_output)
      }
      head(feastdata)
    }


    feastdata1<-feastdata%>%group_by(row)%>%
      summarise_all(sum)
    feastdata1<-feastdata1%>%column_to_rownames(var = "row")
    colsumfeastdata1<-colSums(feastdata1)%>%data.frame()

    feastdata2<-feastdata1
    for (tt in 1:length(colnames(feastdata1))) {
      feastdata2[,tt]<-feastdata1[,tt]/colsumfeastdata1$.[tt]
    }

    groupnum<-match(all_group,rownames(feastdata2))
    groupnum1<-match(rownames(feastdata2),all_group)
    groupnum1[is.na(groupnum1)]<-max(groupnum1[!is.na(groupnum1)])+1
    groupnum[is.na(groupnum)]<-setdiff(groupnum1,groupnum)

    feastdata3<-feastdata2[groupnum,]

    FEAST_output<-feastdata3
  }

  func="times"


  return(list(FEAST_output=FEAST_output,
              sink_group=sink_group,
              source_group=source_group,fun=func,
              otu_table=submchat$otu_table))


}



"calcMicrochatSource" <- function(submchat,
                                  nworker=4,
                                  sink_group = "ml",
                                  source_group="ct") {


  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  sink_group = sink_group
  abun<-submchat$otu_table


  if (sink_group != "all") {

    all_group<-unique(substr(colnames(abun),start = 1,stop = 2))


    stdata<-sourcetracker_file(abun,source_group = source_group,
                               sink_group = sink_group)

    metadata<-stdata$metadata
    otus<-stdata$otus

    data_st<-sourcetrack(metadata,otus,nworker=nworker)

    FEAST_output<-data_st
    colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号

    head(FEAST_output)

  }

  if (sink_group == "all") {

    all_group<-unique(substr(colnames(abun),start = 1,stop = 2))
    samplenum<-length(colnames(abun))

    trid_id<-substr(colnames(abun),start = 1,stop = 2)%>%unique()
    ttk<-sapply(trid_id, function(y){grep(y,colnames(abun))})

    sample.max.size<-lapply(ttk,function(x){length(x)})%>%as.numeric()%>%max()

    feastdata<-data.frame()

    for (tar_group in all_group[2:length(all_group)]) {
      sink_group1=tar_group

      stdata<-sourcetracker_file(abun,source_group = source_group,
                                 sink_group = sink_group1)

      metadata<-stdata$metadata
      otus<-stdata$otus

      data_st<-sourcetrack(metadata,otus,nworker=nworker)

      FEAST_output<-data_st
      colnames(FEAST_output) = paste("repeat_",1:length(colnames(FEAST_output)),sep = "") #取Ids作为平行代号
      FEAST_output<-FEAST_output%>%rownames_to_column(var = "row")

      sample.sie<-length(colnames(FEAST_output))-1
      sample.sie.error<-sample.max.size-sample.sie

      if (sample.sie.error!=0) {
        rowmean.use<-rowMeans(FEAST_output[,2:ncol(FEAST_output)])
        row.na<-length(rowmean.use)
        data.error<-rowmean.use%>%data.frame()
        rownames(data.error)<-rownames(FEAST_output)
        colnames(data.error)<-paste("repeat_",(sample.sie+1):sample.max.size,sep = "")
        FEAST_output<-cbind(FEAST_output,data.error)
      }

      head(FEAST_output)


      {
        feastdata<-rbind(feastdata,FEAST_output)
      }
      head(feastdata)
    }

    feastdata_data1<-rowMeans(feastdata[,2:length(colnames(feastdata))],na.rm = TRUE)%>%data.frame()
    feastdata_data1[,"row"]<-feastdata$row
    colnames(feastdata_data1)[1]<-"value"

    feastdata_data1[,"sink"]<-rep(all_group[2:length(all_group)],each=2)
    FEAST_output<-feastdata_data1
    head(FEAST_output)

  }

  func="common"
  return(list(FEAST_output=FEAST_output,
              sink_group=sink_group,source_group=source_group,fun=func,
              otu_table=submchat$otu_table))
}


"plot_sourcetrack" <- function(microchatSourceobj,
                               color_unknown="red",
                               sourcecolor = colorCustom(5,pal ="gygn" ),
                               plot_type=c("pie","chord","barplot"),
                               repeatiton=TRUE,
                               barplot.nogreyback=FALSE,
                               export_path="sourcetracker") {
  plot_type<-match.arg(plot_type)
  export_path<-paste(export_path,"/data_microbiome/microbial sourcetracking analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  library(dplyr)
  single_st<-microchatSourceobj
  FEAST_output<-single_st$FEAST_output
  sink_group<-single_st$sink_group
  otu_tab<-single_st$otu_table
  all.g<-unique(substr(colnames(otu_tab),start = 1,stop = 2))
  sourcecolor<-sourcecolor[1:length(all.g)]
  names(sourcecolor)<-all.g
  sourcecolor<-c(sourcecolor[match(c(sink_group,single_st$source_group),
                                 names(sourcecolor))],color_unknown)

  ttnum<-length(colnames(FEAST_output))/2
  if (!is.integer(ttnum)) {
    ttnum<-(ttnum%>%as.integer())
  } else {
    ttnum<-ttnum
  }

  if (plot_type=="pie") {
    sourcecolor<-sourcecolor[2:length(sourcecolor)]
    if (sink_group=="all") {
      stop("\n","No pie plot could be seen. Please refer to barplot.")
    }
    if (repeatiton) {
      pdf(paste(export_path,"/sink(",sink_group,")_",
                single_st$source_group,"_",plot_type,".pdf",sep = ""),
          width=10, height=10,family = "serif")
      par(mfrow=c(2,ttnum), mar=c(1,1,1,1))
      for (i in 1:length(colnames(FEAST_output))) {


        labs <- paste0(row.names(FEAST_output)," (", round(FEAST_output[,i]/sum(FEAST_output[,i])*100,2), "%)")
        pie(FEAST_output[,i],labels=labs, init.angle=90,col = sourcecolor,
            border="white",main =paste(sink_group,": ",colnames(FEAST_output)[i],sep = "") )


      }
      dev.off()

      par(mfrow=c(2,ttnum), mar=c(1,1,1,1))
      for (i in 1:length(colnames(FEAST_output))) {


        labs <- paste0(row.names(FEAST_output)," (", round(FEAST_output[,i]/sum(FEAST_output[,i])*100,2), "%)")
        pie(FEAST_output[,i],labels=labs, init.angle=90,col = sourcecolor,
            border="white",main =paste(sink_group,": ",colnames(FEAST_output)[i],sep = "") )


      }


    }else{

      asx = as.data.frame(rowMeans(FEAST_output))
      asx  = as.matrix(asx)
      asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100
      labs <- paste0(row.names(asx_norm)," (", round(asx_norm[,1]/sum(asx_norm[,1])*100,2), "%)")
      pdf(paste(export_path,"/sink(",sink_group,")_",
                single_st$source_group,"_",plot_type,".pdf",sep = ""),
          width=10, height=10,family = "serif")
      pie(asx_norm[,1],labels=labs, init.angle=90,
          col =  sourcecolor,
          border="white",main = sink_group)
      dev.off()

      pie(asx_norm[,1],labels=labs, init.angle=90,
          col =  sourcecolor,
          border="white",main = sink_group)
    }
  }
  if (plot_type=="chord") {
    if (sink_group=="all") {
      stop("\n","No chord plot could be seen. Please refer to barplot.")
    }
    if (repeatiton) {

      for (kk in 1:ncol(FEAST_output)) {
        for (jj in 1:nrow(FEAST_output)) {
          FEAST_output[jj,kk][which(FEAST_output[jj,kk]==0)]<-0.001
        }

      }
      library(graphics)
      circlize::circos.clear()
      pdf(paste(export_path,"/sink(",sink_group,")_",
                single_st$source_group,"_",plot_type,".pdf",sep = ""),
          width=10, height=10,family = "serif")

      par(mfrow=c(2,ttnum), mar=c(1,1,1,1),family="serif")
      for (i in 1:length(colnames(FEAST_output))) {
        asx_norm<-FEAST_output[,i]
        df.plot<-asx_norm%>%data.frame()%>%rownames_to_column(var = "row")
        colnames(df.plot)<-c("row","value")
        df.plot[,"sink"]<-sink_group
        df.plot<-subset(df.plot,select=c(1,3,2))
        df.plot$row<-rownames(FEAST_output)

        group_selnum<-length(unique(df.plot$sink))+length(unique(df.plot$row))
        color.use<-sourcecolor
        grid.color<- color.use[1:group_selnum]
        order.sector<- c(unique(df.plot$sink),unique(df.plot$row))
        names(grid.color)<-order.sector

        circlize::chordDiagram(df.plot,
                               order = order.sector,
                               grid.col = grid.color,
                               link.arr.type="big.arrow",
                               directional=1,
                               small.gap = 1,
                               big.gap = 10,
                               diffHeight=0.03,
                               transparency=0.4,
                               link.visible = TRUE,
                               scale = FALSE,
                               direction.type = "diffHeight+arrows",
                               annotationTrack = c("grid"),
                               reduce = -1,
                               annotationTrackHeight=0.06,
                               preAllocateTracks = list(track.height = 0.1))

        textsize=1;
        lgdxpoi=1;
        lgdypoi=1;
        lgd_title="lgd"

        circlize::circos.track(track.index = 2, panel.fun = function(x, y) {
          xlim = circlize::get.cell.meta.data("xlim")
          xplot = circlize::get.cell.meta.data("xplot")
          ylim = circlize::get.cell.meta.data("ylim")
          sector.name = circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(mean(xlim), ylim[1], sector.name,
                                facing = "clockwise",family = par(family="serif"),
                                niceFacing = TRUE, adj = c(-1, 0.5), cex = textsize)
        }, bg.border = NA)

        lgd <- ComplexHeatmap::Legend(at = names(grid.color),
                                      type = "grid",ncol=1,
                                      legend_gp = grid::gpar(fill = grid.color),
                                      title = lgd_title)
        ComplexHeatmap::draw(lgd, x = unit(0.9, "npc") ,
                             y = unit(0.2, "npc"))

      }
      dev.off()

      circlize::circos.clear()
      par(mfrow=c(2,ttnum), mar=c(1,1,1,1),family="serif")
      for (i in 1:length(colnames(FEAST_output))) {
        asx_norm<-FEAST_output[,i]
        df.plot<-asx_norm%>%data.frame()%>%rownames_to_column(var = "row")
        colnames(df.plot)<-c("row","value")
        df.plot[,"sink"]<-sink_group
        df.plot<-subset(df.plot,select=c(1,3,2))
        df.plot$row<-rownames(FEAST_output)

        group_selnum<-length(unique(df.plot$sink))+length(unique(df.plot$row))
        color.use<-sourcecolor
        grid.color<- color.use[1:group_selnum]
        order.sector<- c(unique(df.plot$sink),unique(df.plot$row))
        names(grid.color)<-order.sector

        circlize::chordDiagram(df.plot,
                               order = order.sector,
                               grid.col = grid.color,
                               link.arr.type="big.arrow",
                               directional=1,
                               small.gap = 1,
                               big.gap = 10,
                               diffHeight=0.03,
                               transparency=0.4,
                               link.visible = TRUE,
                               scale = FALSE,
                               direction.type = "diffHeight+arrows",
                               annotationTrack = c("grid"),
                               reduce = -1,
                               annotationTrackHeight=0.06,
                               preAllocateTracks = list(track.height = 0.1))

        textsize=1;
        lgdxpoi=1;
        lgdypoi=1;
        lgd_title="lgd"

        circlize::circos.track(track.index = 2, panel.fun = function(x, y) {
          xlim = circlize::get.cell.meta.data("xlim")
          xplot = circlize::get.cell.meta.data("xplot")
          ylim = circlize::get.cell.meta.data("ylim")
          sector.name = circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(mean(xlim), ylim[1], sector.name,
                                facing = "clockwise",family = par(family="serif"),
                                niceFacing = TRUE, adj = c(-1, 0.5), cex = textsize)
        }, bg.border = NA)

        lgd <- ComplexHeatmap::Legend(at = names(grid.color),
                                      type = "grid",ncol=1,
                                      legend_gp = grid::gpar(fill = grid.color),
                                      title = lgd_title)
        ComplexHeatmap::draw(lgd, x = unit(0.9, "npc") ,
                             y = unit(0.2, "npc"))

      }

    }else{
      library(graphics)
      circlize::circos.clear()
      pdf(paste(export_path,"/sink(",sink_group,")_",
                single_st$source_group,"_",plot_type,".pdf",sep = ""),
          width=10, height=10,family = "serif")
      par(family="serif")
      asx = as.data.frame(rowMeans(FEAST_output))
      asx  = as.matrix(asx)
      asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100

      df.plot<-asx_norm%>%data.frame()%>%rownames_to_column(var = "row")
      colnames(df.plot)<-c("row","value")
      df.plot[,"sink"]<-sink_group
      df.plot<-subset(df.plot,select=c(1,3,2))

      group_selnum<-length(unique(df.plot$sink))+length(unique(df.plot$row))
      color.use<-sourcecolor
      grid.color<- color.use[1:group_selnum]
      order.sector<- c(unique(df.plot$sink),unique(df.plot$row))
      names(grid.color)<-order.sector

      circlize::chordDiagram(df.plot,
                             order = order.sector,
                             grid.col = grid.color,
                             link.arr.type="big.arrow",
                             directional=1,
                             small.gap = 1,
                             big.gap = 10,
                             diffHeight=0.03,
                             transparency=0.4,
                             link.visible = TRUE,
                             scale = FALSE,
                             direction.type = "diffHeight+arrows",
                             annotationTrack = c("grid"),
                             reduce = -1,
                             annotationTrackHeight=0.06,
                             preAllocateTracks = list(track.height = 0.1))

      textsize=1;
      lgdxpoi=1;
      lgdypoi=1;
      lgd_title="lgd"

      circlize:: circos.track(track.index = 2, panel.fun = function(x, y) {
        xlim = circlize::get.cell.meta.data("xlim")
        xplot = circlize::get.cell.meta.data("xplot")
        ylim = circlize::get.cell.meta.data("ylim")
        sector.name = circlize::get.cell.meta.data("sector.index")
        circlize::circos.text(mean(xlim), ylim[1], sector.name,
                              facing = "clockwise",family = par(family="serif"),
                              niceFacing = TRUE, adj = c(-1, 0.5), cex = textsize)
      }, bg.border = NA)

      lgd <- ComplexHeatmap::Legend(at = names(grid.color),
                                    type = "grid",ncol=1,
                                    legend_gp = grid::gpar(fill = grid.color),
                                    title = lgd_title)
      ComplexHeatmap::draw(lgd, x = unit(0.9, "npc") ,
                           y = unit(0.2, "npc"))
      dev.off()

      circlize::circos.clear()
      par(family="serif")
      asx = as.data.frame(rowMeans(FEAST_output))
      asx  = as.matrix(asx)
      asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100

      df.plot<-asx_norm%>%data.frame()%>%rownames_to_column(var = "row")
      colnames(df.plot)<-c("row","value")
      df.plot[,"sink"]<-sink_group
      df.plot<-subset(df.plot,select=c(1,3,2))

      group_selnum<-length(unique(df.plot$sink))+length(unique(df.plot$row))
      color.use<-sourcecolor
      grid.color<- color.use[1:group_selnum]
      order.sector<- c(unique(df.plot$sink),unique(df.plot$row))
      names(grid.color)<-order.sector

      circlize::chordDiagram(df.plot,
                             order = order.sector,
                             grid.col = grid.color,
                             link.arr.type="big.arrow",
                             directional=1,
                             small.gap = 1,
                             big.gap = 10,
                             diffHeight=0.03,
                             transparency=0.4,
                             link.visible = TRUE,
                             scale = FALSE,
                             direction.type = "diffHeight+arrows",
                             annotationTrack = c("grid"),
                             reduce = -1,
                             annotationTrackHeight=0.06,
                             preAllocateTracks = list(track.height = 0.1))

      textsize=1;
      lgdxpoi=1;
      lgdypoi=1;
      lgd_title="lgd"

      circlize::circos.track(track.index = 2, panel.fun = function(x, y) {
        xlim = circlize::get.cell.meta.data("xlim")
        xplot = circlize::get.cell.meta.data("xplot")
        ylim = circlize::get.cell.meta.data("ylim")
        sector.name = circlize::get.cell.meta.data("sector.index")
        circlize::circos.text(mean(xlim), ylim[1], sector.name,
                              facing = "clockwise",family = par(family="serif"),
                              niceFacing = TRUE, adj = c(-1, 0.5), cex = textsize)
      }, bg.border = NA)

      lgd <- ComplexHeatmap::Legend(at = names(grid.color),
                                    type = "grid",ncol=1,
                                    legend_gp = grid::gpar(fill = grid.color),
                                    title = lgd_title)
      ComplexHeatmap::draw(lgd, x = unit(0.9, "npc") ,
                           y = unit(0.2, "npc"))
    }
  }
  if (sink_group!="all") {
    if (plot_type=="barplot") {
      if (repeatiton) {
        cat("\n","No barplot for this. Please refer to pie plot.","\n")
        message("\n","No barplot for this. Please refer to pie plot.")
      } else{

        feastdata_bar<-FEAST_output%>%rowMeans()%>%data.frame()
        rownames(feastdata_bar)<-rownames(FEAST_output)
        colnames(feastdata_bar)<-"proport"
        feastdata_bar<-feastdata_bar%>%rownames_to_column(var = "row")

        feastdata_bar$label<-paste(round(feastdata_bar$proport*100,2),"%")
        feastdata_bar$proport<-feastdata_bar$proport*100
        feastdata_bar$label.pos<-feastdata_bar$proport+0.03

        feastdata_bar$row<-factor(feastdata_bar$row,levels = unique(feastdata_bar$row))
        names(sourcecolor)[2:length(sourcecolor)]<-c(single_st$source_group,"unknown")
        sourcecolorx<-sourcecolor[match(c(single_st$source_group,"unknown"),
                                       names(sourcecolor))]

        p<-ggplot(feastdata_bar,aes(x=row,y=proport)) +
          geom_col(aes(fill=row),width = 0.5)+
          scale_fill_manual(values = sourcecolorx)+
          scale_y_continuous(expand = c(0,0),limits=c(0,100))+
          geom_text(aes(x = row,y = label.pos,label = label),
                    vjust=-0.5,
                    size = 5,color = "black",family = "serif")+
          labs(x = "Source",y="Proportion (%)",title = paste("Sink ",sink_group,sep = ""))

        if (barplot.nogreyback) {
          p<-p+
            theme_test()+
            theme(#panel.border = element_blank(),
              axis.ticks.length = unit(0.2,"lines"),
              axis.ticks = element_line(color='black'),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = -1),
              axis.title.y=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.text.x=element_text(colour = "black",size = 10,
                                       angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
              strip.text = element_text(colour = "red",size = 10,face = "bold",family = "serif"),
              strip.background = element_rect(fill = "white"),
              #strip.background =  element_blank(),
              legend.position = "none",aspect.ratio = 1)}

        p<-p+theme(#panel.border = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black'),
          #axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = -1),
          axis.title.y=element_text(colour='black', size=12,face = "bold",family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=10,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 10,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          strip.text = element_text(colour = "red",size = 10,face = "bold",family = "serif"),
          strip.background = element_rect(fill = "white"),
          #strip.background =  element_blank(),
          legend.position = "none",aspect.ratio = 1)

        ggsave(paste(export_path,"/sink(",sink_group,")_",single_st$source_group,"_",plot_type,".pdf",sep = ""),p)
        cat("Plot has been exported to ","/",export_path,"",sep = "")
        return(p)
      }
    }
  }
  if (sink_group =="all") {
    if (plot_type=="barplot") {
      if (repeatiton) {
        cat("\n","No barplot for this. Please refer to pie plot.","\n")
        message("\n","No barplot for this. Please refer to pie plot.")
      } else{

        feastdata_bar<-FEAST_output
        feastdata_bar$value.pos<-feastdata_bar$value/2

        for (tt in 1:length(feastdata_bar$row)) {
          if(tt%%2==0) {
            feastdata_bar$value.pos[tt]<-feastdata_bar$value[tt]/2
          }
          if(tt%%2==1) {
            feastdata_bar$value.pos[tt]<-feastdata_bar$value[tt]/2+feastdata_bar$value[tt+1]
          }
        }

        feastdata_bar$label<-paste(round(feastdata_bar$value*100,2),"%")
        feastdata_bar$value<-feastdata_bar$value*100
        feastdata_bar$sink<-factor(feastdata_bar$sink,levels = unique(feastdata_bar$sink))
        names(sourcecolor)[2:length(sourcecolor)]<-c(single_st$source_group,"unknown")
        sourcecolorx<-sourcecolor[match(c(single_st$source_group,"unknown"),
                                        names(sourcecolor))]

        p<- ggplot(feastdata_bar, aes(x=sink, y=value, fill=row)) +
          geom_bar(stat = "identity", width=0.5, col='white') +
          scale_fill_manual(values = sourcecolorx)+
          scale_y_continuous(expand = c(0,0),limits=c(0,100)) +
          geom_text(aes(x = sink,y = value.pos,label = label), vjust=-0.5,
                    size = 5,color = "black",family = "serif")

        if (barplot.nogreyback) {
          p<-p+
            theme_test()+
            theme(#panel.border = element_blank(),
              axis.ticks.length = unit(0.2,"lines"),
              axis.ticks = element_line(color='black'),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=24,face = "bold",family = "serif",vjust = -1),
              axis.title.y=element_text(colour='black', size=24,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_text(colour='black',size=12,family = "serif"),
              axis.text.x=element_text(colour = "black",size = 12,
                                       angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
              strip.text = element_text(colour = "red",size = 10,face = "bold",family = "serif"),
              strip.background = element_rect(fill = "white"),
              #strip.background =  element_blank(),
              legend.position = "none",aspect.ratio = 1)}

        p<-p+ theme(#panel.border = element_blank(),
          text = element_text(family = "serif"),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black'),
          #axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=24,face = "bold",family = "serif",vjust = -1),
          axis.title.y=element_text(colour='black', size=24,face = "bold",family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=12,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 12,
                                   angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
          strip.text = element_text(colour = "red",size = 10,face = "bold",family = "serif"),
          strip.background = element_rect(fill = "white"),
          #strip.background =  element_blank(),
          legend.position = "none",aspect.ratio = 1)

        ggsave(paste(export_path,"/sink(",sink_group,")_",single_st$source_group,"_",plot_type,".pdf",sep = ""),p)
        cat("Plot has been exported to ","/",export_path,"",sep = "")
        return(p)
      }
    }
  }
}


