"randomcolor" <- function(num,select_num=num) {
  r<-round(runif(num, min=0, max=255),0)
  g<-round(runif(num, min=0, max=255),0)
  b<-round(runif(num, min=0, max=255),0)
  randomcolor<-rgb0(r,g,b)
  randomcolor1<-paste(randomcolor[1:select_num],sep = ",")
  return(randomcolor1)
}

"rgb0"<-function(r,g,b){
  return(rgb(r/255,g/255,b/255))
}

"otu_filter" <- function(split_otu,filter_num=6, filter=FALSE,self=FALSE) {

  if (filter) otu <- lapply(split_otu,function(x){

    colname<-colnames(x)

    colname<-substr(colname,start = 1,stop = 2)
    samplenum<-length(colname)
    if (as.integer(samplenum/2)==as.numeric(samplenum/2)) num<-samplenum/2
    if (as.integer(samplenum/2)!=as.numeric(samplenum/2)) num<-as.integer(samplenum/2)-1
    data <- x
    if (!self) num<-num
    if (self) num<-filter_num
    data[data>1] <- 1
    otu <- x[which(rowSums(data) >= num),]
    return(otu)
  })


  if (!filter) otu<-split_otu

  if (!filter) #cat("\n","Samples have been filtered","\n")
  if (filter) #cat("\n","Samples have been filtered","\n")
  return(otu)
}

"data_prep" <- function(mchat) {

  otudata<-mchat$otu_table
  taxdata<-mchat$taxon_table

  taxdata$name<-rownames(taxdata)
  otudata$name<-rownames(otudata)
  tax<-merge(taxdata,otudata,by="name")
  tax<-tax[,1:8]
  rownames(tax)<-tax$name
  tax<-subset(tax,select = -name)
  otu<-subset(otudata,select = -name)

  sampledata<-group_generate(otu)
  rownames(sampledata)<-sampledata$sample
  sampledata<-subset(sampledata,select = -sample)

  micocomobj <- list(otu,tax,sampledata)

  class(micocomobj) <- "microchat"

  names(micocomobj)[[1]]<-"otu_table"
  names(micocomobj)[[2]]<-"taxon_table"
  names(micocomobj)[[3]]<-"sampledata"

  return(list(data=micocomobj,group=unique(sampledata$group)))
}

"otu_split1" <- function(otudata,num=sample_charnum) {

  colname<-colnames(otudata)
  colname<-substr(colname,start = 1,stop = num)
  trt_id <-unique(colname)

  ttk<-sapply(trt_id,function(x){
    grep(x,colnames(otudata))
  })

  split_otu <- lapply(
    apply(ttk,2,FUN = function(x){
      otudata[,x]
    }),function(x){
      if (length(which(rowSums(x)==0))!=0) x[-(which(rowSums(x)==0)),] else x
    })
  return(split_otu)
}

"otu_split2" <- function(otudata,num=sample_charnum) {

  colname<-colnames(otudata)

  colname<-substr(colname,start = 1,stop = num)
  trt_id <-unique(colname)


  ttk<-sapply(trt_id,function(x){
    grep(x,colnames(otudata))
  })

  split_otu <- lapply(
    lapply(ttk,FUN = function(x){
      otudata[,x]
    }),function(x){
      if (length(which(rowSums(x)==0))!=0)  x[-(which(rowSums(x)==0)),] else x
    })
  return(split_otu)
}

"split_merge" <- function(genus_top) {
  otu_rare<-genus_top

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
        rowMeans(x)
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
        rowMeans(x)
      })
  }

  split_otu2<-split_otu%>%data.frame()
  split_otu1<-split_otu%>%data.frame()
  split_otu1<-rownames_to_column(split_otu1,var = "tax")
  split_otu1<-reshape2::melt(split_otu1)
  return(list(split_otu=split_otu2,split_otu1=split_otu1))
}

"group_generate" <- function(otudata,sample_charnum=2) {

  groupname<-substr(colnames(otudata),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
                            grep(x,colnames(otudata))
                                              }),function(y){
                            y %>%length()
                            })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()

if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otu <- otu_split1(otudata,num=sample_charnum)
  } else {
    split_otu <- otu_split2(otudata,num=sample_charnum)
  }

  otu_data<-otu_filter(split_otu)
  treatnum<-length(otu_data)
  names(split_otu)
  samplenum<-length(colnames(otudata))
  sampledata<-data.frame(colnames(otudata),
                         substr(colnames(otudata),start = 1,stop = sample_charnum))
  colnames(sampledata)<-c("sample","group")
  return(sampledata)
}

"random2color" <- function(num,select_num=num,alpha_thres=100) {
  r<-round(runif(num, min=0, max=255),0)
  g<-round(runif(num, min=0, max=255),0)
  b<-round(runif(num, min=0, max=255),0)
  a<-round(runif(num, min=0, max=alpha_thres),0)
  randomcolor<-rgb2(r,g,b,a)
  randomcolor1<-paste(randomcolor[1:select_num],sep = ",")
  return(randomcolor1)
}

"rgb2"<-function(r,g,b,a){
  return(rgb(r/255,g/255,b/255,a/255))
}

"normalize" <- function(data) {
  x<-(data-min(data))/(max(data)-min(data)) ###0-1 standardlization
  y<-(data-mean(data))/(max(data)-min(data))
  z<-(data-mean(data))/sd(data)
  return(list(x=x,y=y,z=z))
}
"scale4dt"<-function(data) {
  mat.use<-as.matrix(data)
  for (ttk in 1:nrow(mat.use)) {
    mat.use[ttk,]<-normalize(mat.use[ttk,])$x
  }
  return(mat.use)
}

"otu_ultra_filter" <- function(abun,filter_num,seed) {
  otu_rare<-abun
  colname<-colnames(otu_rare)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)


  groupname<-substr(colnames(otu_rare),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otu_rare))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()

  if (as.numeric(matchnum) != as.numeric(allnum)) split_otu <- lapply(
    sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),
    FUN = function(x){otu_rare[,x]})

  if (as.numeric(matchnum) == as.numeric(allnum)) split_otu <- apply(
    sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),2,
    FUN = function(x){otu_rare[,x]})


  #if (as.numeric(matchnum) != as.numeric(allnum)) split_otu <- lapply(lapply(sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),FUN = function(x){otu_rare[,x]}),function(x){x[-(which(rowSums(x)==0)),]})

  #if (as.numeric(matchnum) == as.numeric(allnum)) split_otu <- lapply(apply(sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),2,FUN = function(x){otu_rare[,x]}),function(x){x[-(which(rowSums(x)==0)),]})
  set.seed(seed)
  otu<-otu_filter(split_otu,filter_num = filter_num,filter = TRUE,self = TRUE)



  otu <- lapply(otu,function(x){
    x$name <- rownames(x)
    return(x)
  })

  asd<-full_join(otu[[1]],otu[[2]],by=c("name"="name"))

  if (length(otu)>=3) {
    for (tk in 3:length(otu)) {
      asd<-full_join(asd,otu[[tk]],by=c("name"="name"))
    }
  }

  asd<-column_to_rownames(asd,var = "name")
  asd[is.na(asd)]<-0
  return(asd)
}

"colorCustom"  <- function(num,pal=c("pkgn","ywbu","gygn",
                                     "classic","favor",
                                     "cyto_pie","xiena","nicely",
                                     "na1","na2","na3","na4",
                                     "na5","na6","na7","na8",
                                     "set1","set2","set3")) {

  pal<-match.arg(pal)

  if (pal=="nicely") colors<-c("#a58b51" ,"#026250", "#1ea242","#fd2c89","#f9d538",
                               "#09a2fb" ,"#eb84aa","#a6144d", "#d4782b","#aebd28")
  if (pal=="na1") colors<-c("#F4F1DE","#DF7A5E","#3C405B","#82B29A","#F2CC8E")
  if (pal=="na2") colors<-c("#264653","#2A9D8E","#E9C46B","#F3A261","#E66F51")
  if (pal=="na6") colors<-c("#0E606B","#FFC24B","#F66F69","#1597A5","#FEB3AE","#FFF4F2")
  if (pal=="na7") colors<-c("#90C9E7","#219EBC","#136783","#02304A","#FEB705","#FF9E02","#FA8600")

  if (pal=="na3") colors<-c("#B7B5A0","#44757A","#452A3D","#D44C3C","#DD6C4C","#E5855D","#EED5B7")
  if (pal=="na4") colors<-c("#E73847","#F0FAEF","#A8DADB","#457B9D","#1D3557")
  if (pal=="na5") colors<-c("#FBF0C3","#54686F","#E57B7F","#9E3150","#87BBA4")
  if (pal=="na8") colors<-c("#780001","#C11221","#FEF0D5","#002F49","#669BBB")


    if (pal=="favor") colors<-c("#8f850b","#faf4b4","#bde7ff","#dcd9c3","#e9f5fc")
  if (pal=="xiena") colors<-c("#FFC300", "#FF878A" ,"#005C54" ,"#F3F1E4","#AEBD28")
  if (pal=="gygn") colors<-c("#FEBD05", "#889A23" ,"#A69BBA" ,"#0F9DDD", "#887FFE", "#4DB491", "#D4782B",
                             "#3786D4", "#FFC4AB", "#FADD7E", "#8FDDEA", "#83C875", "#DFAB1A", "#31963E",
                             "#659A39" ,"#8EAC0C" ,"#786C06", "#58BEB2", "#98231D" ,"#3AB93E", "#F1D92E",
                             "#840F4F" ,"#9BA4C7", "#369511", "#B00483", "#E22209", "#0B7FE0", "#5BC9EB",
                             "#2D5EF5" ,"#DD8C86", "#38C76E", "#E1B57F" ,"#9DA6EE" ,"#4B4F97", "#C54D80",
                             "#42E95A", "#EAEF76", "#977BC6", "#2674BB" ,"#0E7B72" ,"#9F7423" ,"#D14850",
                             "#A82787", "#C6B342", "#5F015B", "#57FCC7" ,"#69CC7D" ,"#453E52" ,"#CF3949",
                             "#621A0F")

  if (pal=="ywbu") colors<-c("#1BC45E", "#53ACC4", "#FE7808", "#ABB52B", "#605FC3", "#5C4809" ,"#D2A8FE",
                             "#151EF9", "#03DDE6", "#77C30A", "#041295" ,"#D3A81B", "#E7598B", "#0A5748",
                             "#D4F2DE", "#E6BEF5" ,"#D50E68", "#9C87B0", "#68D7E4", "#A05229", "#8E9869",
                             "#F996FE", "#CD0109", "#25F072", "#905B9A", "#EF6E03", "#641688", "#5B0C90",
                             "#E12C2B" ,"#6C1BAD" ,"#70074D", "#8AB98A", "#E7FA56", "#D40435", "#BCCDE4",
                             "#D3E1D2", "#449477", "#A1D5DE", "#A2E5FE", "#48CB53", "#046401", "#51418C",
                             "#A591A6", "#2B74AA" ,"#98D6A1", "#1C791B", "#153906", "#FCD968" ,"#D40D35",
                             "#A658E9")

  if (pal=="pkgn") colors<-c("#1869E9", "#677AB4", "#D86C5A", "#395D44", "#297074" ,"#73946D", "#660D60",
                             "#178AF2", "#EF436E", "#7EA45C", "#4A71A4", "#08319C", "#809606", "#B04E35",
                             "#AC43F7", "#54D41E" ,"#52CBC5", "#D1E12B", "#CFAE1A", "#81F8C0", "#F23E3B",
                             "#6E6B11", "#CA9F0A", "#FA5D1C", "#F0416B", "#C4A3FE", "#089C64", "#241E82",
                             "#0DB9A4" ,"#6D8648", "#00D8A8", "#2CFF29", "#4A8264", "#D6A1E3", "#9C8871",
                             "#50F896", "#6199DA", "#503FA9", "#6AF9B0", "#523D6F" ,"#4C765A", "#7C494B",
                             "#9C9F9C", "#A4C61D", "#9A417C" ,"#A8CEE2", "#C23425", "#4FC4BE", "#45318E",
                             "#15C24A")



  if (pal=="set3") colors<-c("#026250", "#A0FE57", "#56CD70" ,"#0A4CEB" ,"#1EA242",
                             "#F9D538", "#A58B51" ,"#F35DAA" ,"#87B1C9" ,"#A6144D" ,"#EBD408" ,"#561231" ,"#FB2EA8" ,"#7200BE","#EB84AA",
                             "#0ED46D", "#09A2FB", "#1FA36D", "#021510" ,"#164561", "#B02BDC", "#EA7529",
                             "#422321", "#99DE5D" ,"#24683D", "#F8E0D1", "#2DA662" ,"#4CE605" ,"#D6B072",
                             "#0970EB" ,"#D1BD87", "#7D2BD2" ,"#6F0B24" ,"#6601DE")



  if (pal=="set1") colors<-c("#FD2C89" ,"#7AC619", "#026250", "#A0FE57", "#56CD70" ,"#0A4CEB" ,"#1EA242",
                             "#F9D538", "#A58B51" ,"#F35DAA" ,"#87B1C9" ,"#3B21EE" ,"#281C08" ,"#C7BEA1",
                             "#B6100E" ,"#AFB587", "#D266BF" ,"#3A322A", "#C0D0AB" ,"#09373E" ,"#4ACDB2",
                             "#F3EBE2" ,"#CB7391", "#C0AD99" ,"#4E143D" ,"#98FBBA" ,"#76EF3C" ,"#3F57B2",
                             "#C1EBDE", "#8A01AB", "#24683D", "#F8E0D1", "#2DA662" ,"#4CE605" ,"#D6B072",
                             "#0970EB" ,"#D1BD87", "#7D2BD2" ,"#6F0B24" ,"#6601DE")


  if (pal=="set2") colors<-c("#0C8909", "#A6144D" ,"#EBD408" ,"#561231" ,"#FB2EA8" ,"#7200BE","#EB84AA",
                             "#0ED46D", "#09A2FB", "#1FA36D", "#021510" ,"#164561", "#B02BDC", "#EA7529",
                             "#422321", "#99DE5D" ,"#EC359B" ,"#06A566" ,"#298924" ,"#E3E01F" ,"#47AC04",
                             "#C3DF12" ,"#CEC499" ,"#8A7FB2" ,"#189836" ,"#5C1EB1", "#FF927D", "#A30273",
                             "#9436A6" ,"#C90294" ,"#93E698" ,"#D9130A" ,"#24DA32" ,"#A54903" ,"#9F6BC3",
                             "#4D90F2" ,"#DEFCF0" ,"#2090E9", "#1747E6" ,"#EF315E")


  if (pal=="cyto_pie") colors <- c( '#FF7F00', '#984EA3', '#FFED6F','#FFFF33','#377EB8',
                                    '#00D098', '#CCEBC5', '#4DAF4A','#E41A1C','#FDB462',
                                    '#B3DE69', '#BC80BD', '#40ffd1','#B0FF39','#FCCDE5',
                                    '#B02664', '#31a2ff','#ff0000',
                                    '#265BC4', '#513210', '#8300BA','#385E57','#BC226F',
                                    '#657215', '#E04F28', '#2BEFA9','#D9D1ED','#E08294',
                                    '#686632', '#B20077', '#649B7B','#1200C4','#660000',
                                    '#ACBF64', '#254944', '#A36EBA','#79C99B', '#BF738E', '#D8C782','#B72200',
                                    '#FF7F00', '#984EA3', '#FFED6F','#FFFF33','#377EB8',
                                    '#00D098', '#CCEBC5', '#4DAF4A','#E41A1C','#FDB462',
                                    '#B3DE69', '#BC80BD', '#40ffd1','#B0FF39','#FCCDE5',
                                    '#B02664', '#31a2ff','#ff0000',
                                    '#265BC4', '#513210', '#8300BA','#385E57','#BC226F',
                                    '#657215', '#E04F28', '#2BEFA9','#D9D1ED','#E08294',
                                    '#686632', '#B20077', '#649B7B','#1200C4','#660000',
                                    '#ACBF64', '#254944', '#A36EBA','#79C99B', '#BF738E', '#D8C782','#B72200')

  if (pal=="classic") colors<-c("black","white","red","grey","blue")

if (length(colors)<num) {
  colors<-c(colors,colors,colors,colors)
  colors<-colors[1:num]
} else {
  colors<-colors[1:num]
}

  return(colors)
}



"balanceAbun" <- function(otu_table2,maxnum=1.05,minnum=0.95) {
  otu_rare<-otu_table2
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

  ccc1<-calcMean(otu_table2,0)
  ccc2<-calcSd(otu_table2,0)

  add.sample.size=max(samplenum)*2
  fff<-addsample(ccc1,ccc2,0,allnum,add.sample.size,maxnum=maxnum,minnum=minnum)
  newfff<-calcAbun(fff)
  addnum<-add.sample.size-collen%>%as.numeric()
  orinum<-collen%>%as.numeric()
  abunx=otu_table2
  if (length(which(addnum==0))>0) {
    abunxxx<-calcAbun(abunx)
    tkf1<-abunxxx[[1]]
  } else {
    tkf1<-newfff[[1]][,1:addnum[1]]
    colnames(tkf1)<-paste(unique(substr(colnames(tkf1),start = 1,stop = 2)),
                          (1):(addnum[1]),
                          sep="")
  }
  return(tkf1)
}

'adjustTaxaTrend' <- function(tidymchat,
                              specified.taxa="g__Tyzzerella",
                              specified.group=NULL,
                              selected.trend=c(1,1,1,1,1)
) {
  taxon_table<-tidymchat$taxon_table
  otu_table<-tidymchat$otu_table

  allchar<-character()
  for (tk in 1:length(colnames(taxon_table))) {
    allcharr<-taxon_table[,tk]%>%unique()
    {
      allchar<-c(allchar,allcharr)
    }
  }

  if (is.null(specified.taxa)) {
    taxon_table<-taxon_table
  } else {
    specified.taxa1<-specified.taxa[1]
    if (substr(specified.taxa1,start = 2,stop = 3)!="__") {
      specified.taxa<-stringr::str_to_title(specified.taxa)
    } else {
      splitdata<-str_split_fixed(specified.taxa,pattern="__",2)
      splitdata[,2]<-stringr::str_to_title(splitdata[,2])
      specified.taxa<-paste(splitdata[,1],"__",splitdata[,2],sep = "")
    }

    fulltaxon<-lapply(
      lapply(
        sapply(as.list(specified.taxa), grep,allchar),
        function (x) {
          x<-x[1]
        }),
      function(y) {
        allchar[y]
      })%>%as.character()

    taxon_table2<-data.frame()
    for (sk in 1:length(colnames(taxon_table))) {
      taxon_table1<-taxon_table[which(taxon_table[,sk] %in% fulltaxon),]
      {
        taxon_table2<-rbind(taxon_table2,taxon_table1)
      }
    }

    taxon_table<-taxon_table2

    cat("\n","Taxa belonging to '",fulltaxon,"' have (has) been modified.",sep=" ")
  }

  otu_table2<-otu_table[which(rownames(otu_table) %in% rownames(taxon_table)),]

  if (is.null(selected.trend)) {
    if (is.null(specified.group)) {
      abunxxx<-calcAbun(otu_table2)
      tkg1<-balanceAbun(abunxxx[[1]],
                        maxnum=1+0.05,
                        minnum=1-0.05)
      for (ttk in 2:length(abunxxx)) {
        tkg<-balanceAbun(abunxxx[[ttk]],
                         maxnum=1+0.05,
                         minnum=1-0.05)
        {
          tkg1<-cbind(tkg1,tkg)
        }
      }
    } else {
      abunxxx<-calcAbun(otu_table2)
      tkg1<-balanceAbun(abunxxx[[specified.group]],
                        maxnum=1+0.05,
                        minnum=1-0.05)
      for (ttk in setdiff(1:length(abunxxx),specified.group)) {
        tkg<-balanceAbun(abunxxx[[ttk]],
                         maxnum=1+0.05,
                         minnum=1-0.05)
        {
          tkg1<-cbind(tkg1,tkg)
        }
      }

      tkg2<-otu_table2
      for (t in 1:ncol(tkg2)) {
        tkg2[,which(colnames(tkg2)==colnames(tkg1)[t])]<-tkg1[,t]
      }
      tkg1<-tkg2
    }
  } else {
    if (is.null(specified.group)) {
      if (length(selected.trend)<length(unique(group_generate(otu_table2)$group))) {
        stop("Please provide more `selected.trend`")
      } else {
        abunxxx<-calcAbun(otu_table2)
        tkg1<-balanceAbun(abunxxx[[1]],
                          maxnum=selected.trend[1]+0.05,
                          minnum=selected.trend[1]-0.05)
        for (ttk in 2:length(abunxxx)) {
          tkg<-balanceAbun(abunxxx[[ttk]],
                           maxnum=selected.trend[ttk]+0.05,
                           minnum=selected.trend[ttk]-0.05)
          {
            tkg1<-cbind(tkg1,tkg)
          }
        }
      }
    } else {
      if (length(selected.trend)<length(unique(group_generate(otu_table2)$group))) {
        stop("Please provide more `selected.trend`")
      } else {
        abunxxx<-calcAbun(otu_table2)
        tkg1<-balanceAbun(abunxxx[[specified.group]],
                          maxnum=selected.trend[specified.group]+0.05,
                          minnum=selected.trend[specified.group]-0.05)
        for (ttk in setdiff(1:length(abunxxx),specified.group)) {
          tkg<-balanceAbun(abunxxx[[ttk]],
                           maxnum=selected.trend[ttk]+0.05,
                           minnum=selected.trend[ttk]-0.05)
          {
            tkg1<-cbind(tkg1,tkg)
          }
        }

        tkg2<-otu_table2
        for (t in 1:ncol(tkg2)) {
          tkg2[,which(colnames(tkg2)==colnames(tkg1)[t])]<-tkg1[,t]
        }
        tkg1<-tkg2
      }
    }
  }

  for (t in 1:nrow(tkg1)) {
    otu_table[which(rownames(otu_table)==rownames(tkg1)[t]),]<-tkg1[t,]
  }
  tidymchat$otu_table<-otu_table
  tidymchat$taxon_table<-tidymchat$taxon_table
  return(tidymchat)
}

"filter_merge_sum" <- function(otudata,taxon_table) {
  data1<-filter_merge(otudata,taxon_table,
                      filter_num = 1,
                      sample_charnum=2)
  abun<-data1$otu_table
  taxon<-data1$tax_table
  sampledata<-group_generate(abun,sample_charnum=2)

  abun_t<-t(abun)%>%data.frame()
  abun_t$sample<-rownames(abun_t)
  abun_new<-merge(abun_t,sampledata,by="sample")
  abun_new<-subset(abun_new,select = -sample)
  abun_new1<-abun_new%>%
    group_by(group)%>%
    summarise_all(sum)

  abun_new1$group<-factor(abun_new1$group, levels = unique(sampledata$group))
  abun_new1x<-abun_new1
  for (t in 1:nrow(abun_new1x)) {
    abun_new1x[t,]<-abun_new1[match(unique(sampledata$group)[t],abun_new1$group),]
  }
  abun_new1<-abun_new1x

  abun_new1<-tibble::column_to_rownames(abun_new1,var = "group")

  abun_new2<-t(abun_new1)%>%data.frame()
  abun_new4<-abun_new2
  abun_new4[abun_new4!=0]<-1
  abun_new2$name<-rownames(abun_new2)
  abun_new4$name<-rownames(abun_new4)
  taxon$name<-rownames(taxon)

  abun_new3<-merge(abun_new2,taxon,by="name")
  abun_new5<-merge(abun_new4,taxon,by="name")
  return(list(data_venn=abun_new1,
              data_with_taxa=abun_new3,
              data_with_2taxa=abun_new5))
}


"filter_merge" <- function(otudata,taxon_table,
                           filter_num = 4,
                           sample_charnum=2) {

  groupname<-substr(colnames(otudata),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otudata))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()

  if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otu <- otu_split1(otudata,num=sample_charnum)
  } else {
    split_otu <- otu_split2(otudata,num=sample_charnum)
  }

  otu_data<-otu_filter(split_otu,filter_num = filter_num,filter = TRUE,self = TRUE)

  treatnum<-length(otu_data)

  otu_kinds<-NULL
  for (k in 1:treatnum) {
    otu_kind<-rownames(otu_data[[k]])%>%data.frame()

    {
      otu_kinds<-rbind(otu_kinds,otu_kind)
    }

  }
  rownames_summary_afterfilter<-unique(otu_kinds)%>%data.frame()
  colnames(rownames_summary_afterfilter)<-"label"
  rownames_num<-length(rownames_summary_afterfilter$label)

  gmat<-lapply(otu_data, function(x){
    x$name<-rownames(x)
    return(x)
  })


  tttt<-gmat[[1]]
  for (t in 2:treatnum) {
    ddd<-gmat[[t]]
    {
      tttt<-merge(tttt,ddd, all = TRUE)
    }
  }
  if  (rownames_num==length(tttt$name)) cat("\n","data matched，OTU number merged :",rownames_num,"\n")
  if  (rownames_num !=length(tttt$name)) cat("\n","data mismatched，OTU number merged :",rownames_num,"\n")

  tttt[is.na(tttt)]<-0
  rownames(tttt)<-tttt$name
  otu_new<-dplyr::right_join(rownames_summary_afterfilter,tttt,by=c('label'='name'))
  rownames(otu_new)<-otu_new$label
  otu_new<-otu_new[,-1]

  otu_new_data<-otu_new
  otu_new_data$name<-rownames(otu_new_data)
  taxon_table$label<-rownames(taxon_table)
  taxon_newtable<-taxon_table[which(taxon_table$label %in% otu_new_data$name),]
  taxon_newtable<-subset(taxon_newtable,select = -label)

  export_object<-list(network_data=otu_data,otu_table=otu_new,tax_table=taxon_newtable)
  return(export_object)
}


"norm.test"<- function(input.data,alpha=0.05){
  sol<- shapiro.test(input.data)
  if(sol$p.value>alpha){
    cat(paste("success: Obeying normal distribution, p.value=",sol$p.value,">",alpha),"\n")
  }else{
    cat(paste("error  : Not obeying normal distribution, p.value=",sol$p.value,"<=",alpha),"\n")
  }
  sol
}

"calcSegmentOrder" <- function(data.alpha,y.ratio,data_poi) {
  data.seg<-data.alpha

  data.seg<-subset(data.seg, select=c(index,group1,group2,p.value,sig,yvalue1,yvalue2,order,error))
  data.seg$order1<-match(data.seg$group1,unique(data_poi$group))
  data.seg$order2<-match(data.seg$group2,unique(data_poi$group))

  data.seg$y1<-data.seg$yvalue1
  data.seg$y2<-data.seg$yvalue2

  data.segment<-data.frame(
    group=c(data.seg$group1,data.seg$group2),
    yvalue=c(data.seg$yvalue1,data.seg$yvalue2),
    y=c(data.seg$y1,data.seg$y2))

  data.segment$order<-rownames(data.segment)
  data.segment.all<-data.frame()

  for (uname in unique(data.segment$group)) {
    data.segment1<- data.segment[which(data.segment$group==uname),]

    for (tt in 1:nrow(data.segment1)) {

      data.segment1$y[tt]<-data.segment1$yvalue[tt]+(tt)*y.ratio*2
      {
        data.segment.all<-rbind(data.segment.all,data.segment1[tt,])
      }
    }
  }

  data.segment.all<-data.segment.all[order(data.segment.all$order,decreasing = FALSE),]
  group.num<-nrow(data.segment.all)
  group.num.half<-group.num/2

  data.seg$y1<-data.segment.all$y[1:group.num.half]
  data.seg$y2<-data.segment.all$y[(group.num.half+1):group.num]
  data.seg<-subset(data.seg,select=c(index,group1,group2,p.value,sig,order1,order2,
                                     y1,y2))

  data.seg$y3<-data.seg$y1+y.ratio
  data.seg$y4<-data.seg$y2+y.ratio


  for (kk in 1:nrow(data.seg)) {
    if (data.seg$y3[kk]<data.seg$y4[kk]) {
      data.seg$y3[kk]<-data.seg$y4[kk]
    } else {
      data.seg$y3[kk]<-data.seg$y3[kk]
      data.seg$y4[kk]<-data.seg$y3[kk]
    }
  }


  data.segment<-data.frame(
    group=c(data.seg$group1,data.seg$group2),
    ystart=c(data.seg$y1,data.seg$y2),
    ystop=c(data.seg$y3,data.seg$y4)
  )

  data.segment$order<-rownames(data.segment)
  data.segment.all<-data.frame()

  for (uname in unique(data.segment$group)) {
    data.segment1<- data.segment[which(data.segment$group==uname),]

    if (nrow(data.segment1)>=2) {
      for (tt in 2:nrow(data.segment1)) {
        if (data.segment1$ystart[tt]<data.segment1$ystop[tt-1]){
          data.segment1$ystart[tt]<-data.segment1$ystop[tt-1]+y.ratio
        }}

    }
    {
      data.segment.all<-rbind(data.segment.all,data.segment1)
    }
  }


  data.segment.all<-data.segment.all[order(data.segment.all$order,decreasing = FALSE),]
  group.num<-nrow(data.segment.all)
  group.num.half<-group.num/2

  data.seg$y1<-data.segment.all$ystart[1:group.num.half]
  data.seg$y2<-data.segment.all$ystart[(group.num.half+1):group.num]
  data.seg$y3<-data.segment.all$ystop[1:group.num.half]
  data.seg$y4<-data.segment.all$ystop[(group.num.half+1):group.num]

  data.seg<-subset(data.seg,select=c(index,group1,group2,p.value,sig,order1,order2,
                                     y1,y2,y3,y4))

  return(data.seg)
}

"calcComparisonOrder" <- function(my_comparisons,data_poi) {


  newdata<-data.frame()
  for (tt in 1:length(my_comparisons)) {
    if (tt==1) {
      my_comparisons1<-my_comparisons[tt]%>%data.frame()
      {
        newdata<-rbind(newdata,my_comparisons1)
      }
    } else {
      my_comparisons1<-my_comparisons[tt]%>%data.frame()

      {
        newdata<-cbind(newdata,my_comparisons1)
      }
    }
  }
  my_comparisons1<-newdata%>%t()%>%data.frame()
  colnames(my_comparisons1)<-c("group1","group2")
  rownames(my_comparisons1)<-paste("comparison",1:nrow(my_comparisons1),sep = "")

  group.selected<-c(my_comparisons1$group1,
                    my_comparisons1$group2)

  muti.index<-group.selected[duplicated(group.selected)]

  for (tt in 1:length(my_comparisons)) {
    new.com<-my_comparisons[[tt]]
    if (is.na(match(my_comparisons1$group1[tt],muti.index)) &
        (my_comparisons1$group2[tt] %in% muti.index)) {
      my_comparisons1$group1[tt]<-new.com[2]
      my_comparisons1$group2[tt]<-new.com[1]
    }
  }

  my_comparisons1
  my_comparisons2<-my_comparisons1%>%t()%>%data.frame()%>%as.list()
  my_comparisons<-my_comparisons2



  sta<-lapply(my_comparisons, function(x){
    alphadiv.use.selected1<-data_poi[
      which(data_poi$group==x[1] | data_poi$group==x[2]),]
    stats<-data.frame(
      index=index,
      group1=x[1],
      group2=x[2],
      max1=max(alphadiv.use.selected1$value[
        which(alphadiv.use.selected1$group==x[1])]),
      max2=max(alphadiv.use.selected1$value[
        which(alphadiv.use.selected1$group==x[2])])
    )
  })

  data.alpha<-data.frame()
  for (tt in 1:length(sta)) {
    data.al<-sta[[tt]]
    {
      data.alpha<-rbind(data.alpha,data.al)
    }
  }

  data.alpha1<-data.alpha[order(data.alpha$max1,data.alpha$max2,decreasing = FALSE),]
  data.alpha2<-subset(data.alpha1, select=c(group1,group2))%>%t()%>%data.frame()%>%as.list()
  return(data.alpha2)
}

"variance.test"<- function(input.data,index,alpha=0.05){
  sol<- bartlett.test(as.formula(paste(index, "~ group",sep = "")), data = input.data)
  if (!is.na(sol$p.value)) { if(sol$p.value>alpha){
    cat(paste("success: Showing homogeneity of variance, p.value=",sol$p.value,">",alpha),"\n")
  }else{
    cat(paste("error  : Not showing homogeneity of variance, p.value=",sol$p.value,"<=",alpha),"\n")
  }
    }
  sol
}

"aov.test"<- function(input.data,index){
  fit<- aov(as.formula(paste(index, "~ group",sep = "")), data = input.data)
}


"kruskal.pre"<- function(input.data,index){
  fit<- kruskal.test(as.formula(paste(index, "~ group",sep = "")), data = input.data)
}

"dunn.test"<- function(input.data,index){
  fit<- rstatix::dunn_test(as.formula(paste(index, "~ group",sep = "")), data = input.data, p.adjust.method = "bonferroni")
}


"mylr" = function(x,y){
  x_mean = mean(x)
  y_mean = mean(y)
  xy_mean = mean(x*y)
  xx_mean = mean(x*x)
  yy_mean = mean(y*y)

  m = (x_mean*y_mean - xy_mean)/(x_mean^2 - xx_mean)
  b = y_mean - m*x_mean

  f = m*x+b# 线性回归方程

  sst = sum((y-y_mean)^2)
  sse = sum((y-f)^2)
  ssr = sum((f-y_mean)^2)

  result = c(m,b,sst,sse,ssr)
  names(result) = c('m','b','sst','sse','ssr')

  return(result)
}



"color_match" <- function(otudata,taxdata,tree11,
                          taxa_order=taxa_order,taxa_num=taxa_num) {
  tree_filter<-tree_filter(tree11=tree11,taxon=taxdata,abun=otudata)
  color_group<-color_mannual_use(tree_filter$data_tree)
  color_group<-color_group%>%as.data.frame()
  color_group$taxa<-rownames(color_group)
  color_group$taxa<-substr(color_group$taxa,start = 4,stop = str_length(color_group$taxa))

  data_prep<-data_prep(otudata,taxdata)
  data<-data_prep$data
  group_order<-data_prep$group

  t1 <- trans_abund$new(dataset = data,
                        taxrank = names(data$taxa_abund)[taxa_order],
                        ntaxa =taxa_num)

  taxa_select<-data.frame(t1$data_taxanames,rep(0,length(t1$data_taxanames)))
  colnames(taxa_select)<-c("taxon","zero")

  color_group_select<-color_group[which(color_group$taxa %in% t1$data_taxanames),]
  color_group1_select<-dplyr::right_join(taxa_select,color_group_select,by=c('taxon'='taxa'))
  colnames(color_group1_select)[3]<-"color"

  return(color_group1_select$color)
}


"filepack" <- function(file) {
  usethis::use_data(file,overwrite = TRUE)
}



"calcSegment2Order" <- function(data.alpha,y.ratio,data_poi) {
  data.seg<-data.alpha

  data.seg<-subset(data.seg, select=c(index,group1,group2,p.value,sig,yvalue1,yvalue2,order,error))
  data.seg$order1<-match(data.seg$group1,unique(data_poi$group))
  data.seg$order2<-match(data.seg$group2,unique(data_poi$group))

  data.seg$y1<-data.seg$yvalue1
  data.seg$y2<-data.seg$yvalue2

  data.seg<-subset(data.seg,select=c(index,group1,group2,p.value,sig,order1,order2,
                                     y1,y2))



  data.seg$order.error<-abs(data.seg$order2-data.seg$order1)
  data.seg$yvalue.error<-abs(data.seg$y2-data.seg$y1)
  data.seg<-data.seg[order(data.seg$order.error,
                           data.seg$yvalue.error),]

  data.seg<-subset(data.seg,select=c(index,group1,group2,p.value,sig,order1,order2,
                                     y1,y2))


  data.seg$y3<-data.seg$y1+y.ratio
  data.seg$y4<-data.seg$y2+y.ratio


  for (kk in 1:nrow(data.seg)) {
    if (data.seg$y3[kk]<data.seg$y4[kk]) {
      data.seg$y3[kk]<-data.seg$y4[kk]
    } else {
      data.seg$y3[kk]<-data.seg$y3[kk]
      data.seg$y4[kk]<-data.seg$y3[kk]
    }
  }

  data.seg$y1<-data.seg$y1+y.ratio
  data.seg$y2<-data.seg$y2+y.ratio
  data.seg$y3<-data.seg$y3+y.ratio
  data.seg$y4<-data.seg$y4+y.ratio



  if (nrow(data.seg)>1) {

    for (tk in 2:nrow(data.seg)) {
      data.seg1<-data.seg[1:tk,]
      sink.char<-data.seg1$group1[tk]


      match.num1<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group1)
      })%>%as.numeric()

      match.num1x<-match.num1%>%max()
      len1<-match.num1%>%length()


      match.num2<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group2)
      })%>%as.numeric()

      match.num2x<-match.num2%>%max()
      len2<-match.num2%>%length()

      if (is.na(match.num1x) & !is.na(match.num2x)) {
        data.seg[tk,]$y1<-data.seg$y2[match.num2x]+y.ratio*2
      }
      if (!is.na(match.num1x) & is.na(match.num2x)) {
        data.seg[tk,]$y1<-data.seg$y1[match.num1x]+y.ratio*2
      }
      if (!is.na(match.num1x) & !is.na(match.num2x)) {
        if (match.num1x<match.num2x) {
          match.num.row<-match.num2x
          data.seg[tk,]$y1<-data.seg$y2[match.num.row]+y.ratio*2
        } else {
          match.num.row<-match.num1x
          data.seg[tk,]$y1<-data.seg$y1[match.num.row]+y.ratio*2
        }
      }



    }


    for (tk in 2:nrow(data.seg)) {
      data.seg1<-data.seg[1:tk,]
      sink.char<-data.seg1$group2[tk]


      match.num1<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group1)
      })%>%as.numeric()

      match.num1x<-match.num1%>%max()
      len1<-match.num1%>%length()


      match.num2<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group2)
      })%>%as.numeric()

      match.num2x<-match.num2%>%max()
      len2<-match.num2%>%length()

      if (is.na(match.num1x) & !is.na(match.num2x)) {
        data.seg[tk,]$y2<-data.seg$y2[match.num2x]+y.ratio*2
      }
      if (!is.na(match.num1x) & is.na(match.num2x)) {
        data.seg[tk,]$y2<-data.seg$y1[match.num1x]+y.ratio*2
      }
      if (!is.na(match.num1x) & !is.na(match.num2x)) {
        if (match.num1x<match.num2x) {
          match.num.row<-match.num2x
          data.seg[tk,]$y2<-data.seg$y2[match.num.row]+y.ratio*2
        } else {
          match.num.row<-match.num1x
          data.seg[tk,]$y2<-data.seg$y1[match.num.row]+y.ratio*2
        }
      }
    }

    for (tk in 2:nrow(data.seg)) {
      data.seg1<-data.seg[1:tk,]
      sink.char<-data.seg1$group1[tk]


      match.num1<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group1)
      })%>%as.numeric()

      match.num1x<-match.num1%>%max()
      len1<-match.num1%>%length()


      match.num2<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group2)
      })%>%as.numeric()

      match.num2x<-match.num2%>%max()
      len2<-match.num2%>%length()

      if (is.na(match.num1x) & !is.na(match.num2x)) {
        data.seg[tk,]$y1<-data.seg$y2[match.num2x]+y.ratio*2
      }
      if (!is.na(match.num1x) & is.na(match.num2x)) {
        data.seg[tk,]$y1<-data.seg$y1[match.num1x]+y.ratio*2
      }
      if (!is.na(match.num1x) & !is.na(match.num2x)) {
        if (match.num1x<match.num2x) {
          match.num.row<-match.num2x
          data.seg[tk,]$y1<-data.seg$y2[match.num.row]+y.ratio*2
        } else {
          match.num.row<-match.num1x
          data.seg[tk,]$y1<-data.seg$y1[match.num.row]+y.ratio*2
        }
      }



    }


    for (tk in 2:nrow(data.seg)) {
      data.seg1<-data.seg[1:tk,]
      sink.char<-data.seg1$group2[tk]


      match.num1<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group1)
      })%>%as.numeric()

      match.num1x<-match.num1%>%max()
      len1<-match.num1%>%length()


      match.num2<-sapply(sink.char,function(x) {
        grep(x,data.seg1[1:(tk-1),]$group2)
      })%>%as.numeric()

      match.num2x<-match.num2%>%max()
      len2<-match.num2%>%length()

      if (is.na(match.num1x) & !is.na(match.num2x)) {
        data.seg[tk,]$y2<-data.seg$y2[match.num2x]+y.ratio*2
      }
      if (!is.na(match.num1x) & is.na(match.num2x)) {
        data.seg[tk,]$y2<-data.seg$y1[match.num1x]+y.ratio*2
      }
      if (!is.na(match.num1x) & !is.na(match.num2x)) {
        if (match.num1x<match.num2x) {
          match.num.row<-match.num2x
          data.seg[tk,]$y2<-data.seg$y2[match.num.row]+y.ratio*2
        } else {
          match.num.row<-match.num1x
          data.seg[tk,]$y2<-data.seg$y1[match.num.row]+y.ratio*2
        }
      }
    }

  } else {
    stop("......")
  }

  data.seg$y3<-data.seg$y1+y.ratio
  data.seg$y4<-data.seg$y2+y.ratio


  for (kk in 1:nrow(data.seg)) {
    if (data.seg$y3[kk]<data.seg$y4[kk]) {
      data.seg$y3[kk]<-data.seg$y4[kk]
    } else {
      data.seg$y3[kk]<-data.seg$y3[kk]
      data.seg$y4[kk]<-data.seg$y3[kk]
    }
  }

  for (ktt in 2:nrow(data.seg)) {
    data.segx<-data.seg[ktt,]
    if (data.segx$y3<data.seg[(ktt-1),]$y3) {
      ori<-c(data.seg[(ktt-1),]$order1:data.seg[(ktt-1),]$order2)
      sinkvec<-c(data.segx$order1:data.segx$order2)

      intevec<-intersect(ori,sinkvec)%>%length()
      if (intevec==1 | intevec==0 ) {

        data.seg[ktt,]$y1<-data.seg[ktt,]$y1
        data.seg[ktt,]$y2<-data.seg[ktt,]$y2

      } else if (intevec==2) {
        data.seg[ktt,]$y1<-data.seg[(ktt-1),]$y3+y.ratio
        data.seg[ktt,]$y2<-data.seg[(ktt-1),]$y3+y.ratio
        data.seg[ktt,]$y3<-max(c(data.seg[ktt,]$y1,data.seg[ktt,]$y2))+y.ratio
      }
    } else {
      if (data.segx$y1<data.seg[(ktt-1),]$y3 ) {
        data.seg[ktt,]$y1<-data.seg[(ktt-1),]$y3+y.ratio
      }
      if (data.segx$y2<data.seg[(ktt-1),]$y3) {
        data.seg[ktt,]$y2<-data.seg[(ktt-1),]$y3+y.ratio
      }

    }

  }





  return(data.seg)
}



"scale1dt" <- function(data) {
  mat.use<-as.matrix(data)
  for (ttk in 1:ncol(mat.use)) { mat.use[,ttk]<-normalize(mat.use[,ttk])$x }
  return(mat.use)
}

"scale2dt" <- function(data) {
  mat.use<-as.matrix(data)
  for (ttk in 1:ncol(mat.use)) { mat.use[,ttk]<-normalize(mat.use[,ttk])$y }
  return(mat.use)
}

"scale3dt" <- function(data) {
  mat.use<-as.matrix(data)
  for (ttk in 1:ncol(mat.use)) { mat.use[,ttk]<-normalize(mat.use[,ttk])$z }
  return(mat.use)
}


"kruskalmuticomp" <- function(data,group,compare,value,
                              reversed = FALSE){
  library(multcompView)
  library(pgirmess)
  library(multcomp)
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)

    k <- kruskalmc(value ~ g1, data=sub_dat, probs=0.05)
    dif <- k$dif.com$stat.signif
    names(dif) <- rownames(k$dif.com)
    difL <- multcompLetters(dif,reversed = reversed)
    label <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    label$compare = rownames(label)
    label$type <- i

    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,label,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}



"microchatParamSink" <- function(FileNamePrefix="param_ss_") {

  filesdir<-paste(FileNamePrefix,"*.txt",sep = "")

  datafiles<-lapply(Sys.glob(filesdir),read.delim)
  myfiles0 <- Sys.glob(filesdir)

  myfilesx<-strsplit(myfiles0,'.txt')%>%as.character()
  snlen<-nchar(FileNamePrefix)
  smlen<-nchar(myfilesx)%>%max()
  newname<-substr(myfilesx,start = snlen+1,stop = smlen)

  for (tk in 1:length(datafiles)) {
    names(datafiles)[tk]<-newname[tk]
    datafiles[[tk]]<-datafiles[[tk]]%>%tibble::column_to_rownames(var = "X")
  }

  return(datafiles)
}



"groupOTU.phylo" <- function(phy, focus, group_name="group", ...) {
  attr(phy, group_name) <- NULL
  if ( is(focus, "list") ) {
    for (i in 1:length(focus)) {
      phy <- microchat::gfocus(phy, focus[[i]], group_name, names(focus)[i], ...)
    }
  } else {
    phy <- microchat::gfocus(phy, focus, group_name, ...)
  }
  res <- attr(phy, group_name)
  res[is.na(res)] <- 0
  attr(phy, group_name) <- factor(res)
  return(phy)
}


"gfocus" <- function(phy, focus, group_name, focus_label=NULL, overlap="overwrite") {
  overlap <- match.arg(overlap, c("origin", "overwrite", "abandon"))

  if (is.character(focus)) {
    focus <- which(phy$tip.label %in% focus)
  }

  n <- treeio::getNodeNum(phy)
  if (is.null(attr(phy, group_name))) {
    foc <- rep(0, n)
  } else {
    foc <- attr(phy, group_name)
  }
  i <- max(suppressWarnings(as.numeric(foc)), na.rm=TRUE) + 1
  if (is.null(focus_label)) {
    focus_label <- i
  }

  ## sn <- phy$edge[which.edge(phy, focus),] %>% as.vector %>% unique
  hit <- unique(as.vector(phy$edge[which.edge(phy, focus),]))
  if (overlap == "origin") {
    sn <- hit[is.na(foc[hit]) | foc[hit] == 0]
  } else if (overlap == "abandon") {
    idx <- !is.na(foc[hit]) & foc[hit] != 0
    foc[hit[idx]] <- NA
    sn <- hit[!idx]
  } else {
    sn <- hit
  }

  if (length(sn) > 0) {
    foc[sn] <- focus_label
  }

  attr(phy, group_name) <- foc
  phy
}





"sim_rma_data" <- function(n, k, means = NULL, poly_order = NULL, noise_sd = 10, between_subject_sd = 40, NAs = 0) {

  # Create data structure and simulate data --------------------------------------------------

  # Create empty n x k matrix
  rma_data = matrix(NA, nrow = n, ncol = k + 1)

  # Add column with subject_id
  rma_data[, 1] = 1:n

  if (!is.null(means)) {
    con_means = means

    # Check if length of mean vector corresponds to k
    if (length(means) != k) {
      k = length(means)
      print("Number of factors (k) was changed, because the length of means vector and argument k do not correspond.")
    }
  } else {

    # Simulate conditional means
    if (is.null(poly_order)) {
      con_means = runif(k, min = 100, max = 300)
    } else {

      # Generate polinomial conditional means
      factors   = runif((poly_order + 1), min = 0, max = 1)
      x         = order(runif(k, min = 100, max = 300), decreasing = FALSE)
      con_means = matrix(factors[1], nrow = k)

      for (p in (2:(poly_order + 1))) {
        con_means = con_means + factors[p] * x^p
      }
    }
  }

  # Add con_mean to the rma_data matrix
  rma_data[, 2:(k + 1)] = matrix(rep(con_means, each = n), nrow = n)

  # Simulate subject means Calculate the deviation from the conditional mean for each subject
  mean_deviation        = rnorm(n, mean = 0, sd = between_subject_sd)
  rma_data[, 2:(k + 1)] = rma_data[, 2:(k + 1)] + mean_deviation

  # Check if only one noise_sd was passed (no sphericity) and create vector
  if (length(noise_sd) == 1) {
    noise_sd = rep(noise_sd, times = k)
  }

  # Error if noise_sd vector not a vector of length k
  if (length(noise_sd) != k) {
    print("The vector passed for noise_sd does not have the length k. Please pass a vector of the length k.")
    return(NULL)
  }


  # Add noise to data ----------------------------------------------------------

  noise = matrix(NA, nrow = n, ncol = k)

  for (i in 1:k) {
    noise[, i] = rnorm(n, mean = 0, sd = noise_sd[i])
  }

  rma_data[, 2:(k + 1)] = rma_data[, 2:(k + 1)] + noise

  # Simulating NAs, adds the number of NA to the the data passed as argument
  if (NAs > 0) {
    for (i in 1:NAs) {
      rma_data[runif(1, min = 1, max = n), runif(1, min = 2, max = (k + 1))] = NA
    }
  }


  # Naming columns ------------------------------------------------------------

  factor_names    = character(k + 1)
  factor_names[1] = "Subject_id"

  for (i in 1:k) {
    factor_names[i + 1] = paste("Factor", i)
  }
  colnames(rma_data) = factor_names

  return(data.frame(rma_data))
}

"rma_opc" <- function(rma_data, id = 1, maxpoly = NA, print_plot = FALSE) {

  # Check if the data meet the following requirement:

  # id must be an integer specifying the column position of the ID variable
  if (id %in% 1:ncol(rma_data) == FALSE || length(id) != 1) {
    stop("id must be an integer specifying the column position of the ID variable")
  }

  dependent_variable = as.matrix(rma_data[, -id])


  # Define some variables ---------------------------------------------------

  # number of entities
  n = nrow(rma_data)

  # number of factor levels
  k = ncol(dependent_variable)

  # Specify the names of the 'id'-variable and of the 'condition'-variables
  rm_names = colnames(dependent_variable)
  id_names = colnames(rma_data)[id]


  # check if the data meet the requirements ---------------------------------

  # rma_data needs to meet the following requirements:

  # all variables must be numeric
  if (all(sapply(rma_data, is.numeric)) == FALSE | any(sapply(rma_data, is.factor))) {
    stop("All variables in rma_data must be numeric")
  }

  # n > k (i.e. more entities than factor levels)
  # if (n <= k) {
  #  message("Number of entities must exceed number of factor levels")
  # }

  # k >= 2 (i.e. at least two or more factor levels)
  if (k < 2) {
    stop("At least two factor factor levels required")
  }


  # Convert to long format --------------------------------------------------

  rma_data_long = reshape(rma_data,
                          varying       = rm_names,
                          v.names       = "value",
                          timevar       = "condition",
                          times         = (1:k),
                          idvar         = id_names,
                          new.row.names = 1:(k * n),
                          direction     = "long")
  colnames(rma_data_long)[1] = "id"
  rma_data_long$condition    = as.numeric(rma_data_long$condition)


  # Define some more variables ----------------------------------------------

  # factor level means
  Flm = tapply(rma_data_long$value, rma_data_long$condition, mean)

  # general mean
  Gm = mean(rma_data_long$value)

  # entity/subject mean
  Em = tapply(rma_data_long$value, rma_data_long$id, mean)

  # Measurements
  Me = 1:k

  # Mean of each measurement condition ('MeFlm' dataframe)
  MeFlm     = data.frame(Me, Flm)
  MeFlmlong = MeFlm[rep(seq_len(nrow(MeFlm)), each = n), ]

  # Entities
  E = 1:n

  # Mean of each entity/subject ('EEm' dataframe)
  EEm     = data.frame(E, Em)
  EEmlong = EEm[rep(seq_len(nrow(EEm)), each = k), ]


  # Orthogonal polynomial Contrasts -----------------------------------------

  # maximal polynomial degree for orthogonal polynomials
  if ((maxpoly > k - 1) | (is.na(maxpoly))) {
    maxpoly = k - 1
  }

  # Defining Contrast weights for orthogonal polynomial contrasts
  contrast_weights = (t(contr.poly(k)))[1:maxpoly, ]

  # Applying formula for linear contrasts
  weighted_dependend_variables = dependent_variable[rep(1:n, each = maxpoly), ] * (contrast_weights)[rep(1:maxpoly, n), ]
  linear_subject_contrasts     = matrix(rowSums(weighted_dependend_variables), byrow = TRUE, ncol = maxpoly)

  # Computing contrast estimators for each orthogonal polynomial contrast as well as standard errors for thees estimators
  contrast_estimator = colMeans(linear_subject_contrasts)
  contrast_se        = sqrt(apply(linear_subject_contrasts, 2, var)) / sqrt(n)

  # Computing t-values for each contrast
  contrast_t_values = contrast_estimator / contrast_se
  # contrast_F_values = contrast_t_values^2

  # Computing the corresponding p-values
  contrast_p_values = 1 - pt(abs(contrast_t_values), n - 1)

  # Computing sums of squares for each contrast
  contrast_ss = n * contrast_estimator^2 / rowSums(contrast_weights^2)

  # Computing amount of the variance in the dependent variable explained by the factor which in turn can be explained by a cerain
  # orthogonal polynomial trend ss_trend / ss_factor
  proportional_trend_contribution = contrast_ss / rep(sum(rep((Flm - Gm)^2, each = n)), maxpoly)


  # Create contrast table ---------------------------------------------------

  # define source variable
  source         = rownames(contrast_weights)
  contrast_table = data.frame(check.names = FALSE,
                              Source      = source,
                              `Sum of squares`     = contrast_ss,
                              `Proportional contribution to the factor effect` = proportional_trend_contribution,
                              `Contrast estimator` = contrast_estimator,
                              `Standard error`     = contrast_se,
                              `Degrees of freedom` = rep((n - 1), maxpoly),
                              `t-value` = contrast_t_values,
                              `p-value` = contrast_p_values)
  rownames(contrast_table) = NULL


  # Orthogonal polynomial trends as polynomial regression ------------------
  # Used to plot the fitted polynomials This is used to display the aditional explanation of the variance in the dependent
  # variable by adding higher order trendcomponent successively

  # Initialize empty dataframe for polynomial regression coefficients
  poly_coef = matrix(0, ncol = maxpoly, nrow = k)

  # Fitting the orthogonal Polynomials In each cycle of the loop the coefficients are assigned to the i-th column of the object
  # poly_coef
  for (i in 1:maxpoly) {
    poly                      = lm(rma_data_long$value ~ poly(rma_data_long$condition, degree = i, raw = TRUE))
    poly_coef[, i][1:(i + 1)] = poly$coef
  }


  # Plotting contrats (ggplot) -------------------------------------------------

  # create datapoints for polynomial plot: this code automatically sets up the data that is required to plot the k-1 polynomial
  # regression lines
  x               = seq(1, k, length.out = 1000)
  poly_values     = outer(x, 0:(k - 1), `^`)
  poly_curve_data = data.frame(x, poly_values %*% poly_coef)
  poly_curve_data = gather(poly_curve_data, curve, y, -x)

  # plot the k-1 polynomial regression lines
  poly_plot = ggplot(data = rma_data_long, aes(x = condition, y = value)) +
    geom_point() +
    labs(col   = "Order of \npolynomial",
         x     = "Condition",
         y     = "Value",
         title = "Orthogonal polynomial contrasts") +
    geom_path(data = poly_curve_data, aes(x, y, color = curve), lwd = 1.2) +
    scale_color_discrete(labels = as.character(1:(maxpoly))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key       = element_rect(colour = "black"),
          plot.title       = element_text(face = "bold", hjust = 0.5))


  # Return the contrast-table and plot ----------------------------------------

  if (print_plot == TRUE) {
    print(poly_plot)
  }

  return(list(contrast_table = contrast_table, poly_plot = poly_plot,rma_data_long=rma_data_long,poly_curve_data=poly_curve_data))
}



"find_overlap" <- function(mods, bigtable){
  vec1 = bigtable[,which(names(bigtable)==mods[1])]
  vec2 = bigtable[,which(names(bigtable)==mods[2])]
  return(sum(vec1*vec2==1))
}

"find_only_in_1" <- function(mods, bigtable){
  vec1 = bigtable[,which(names(bigtable)==mods[1])]
  vec2 = bigtable[,which(names(bigtable)==mods[2])]
  return(sum(vec1==1 & vec2==0))
}

"find_only_in_2" <- function(mods, bigtable){
  vec1 = bigtable[,which(names(bigtable)==mods[1])]
  vec2 = bigtable[,which(names(bigtable)==mods[2])]
  return(sum(vec2==1 & vec1==0))
}

"find_N" <- function(mods, mapping, bigtable){
  nwk1 = bigtable[, which(mapping$Group== mapping $Group[which(mapping $group==mods[1])])]
  nwk2 = bigtable[, which(mapping$Group== mapping $Group[which(mapping $group==mods[2])])]
  match_nwk1_nwk2 = sum((rowSums(nwk1) + rowSums(nwk2))>0)
  return(match_nwk1_nwk2)
}

"fisher_test" <- function(x){
  contingency_table <- matrix(unlist(matrix(data.frame(x[3:6]), nrow=2)), nrow=2)
  test_p = fisher.test(contingency_table, alternative = "greater")$p.value
  return(test_p)
}

"model_preserve" <- function(
    node_table2 = node_table2,
    n = 3,
    padj = FALSE,
    export_path="microbial network analysis"
){

  require(tidyfst)
  dir.create(export_path, recursive = TRUE)
  node_table2$value = 1
  #提取全部的模块及其名称
  module_list = unique(node_table2$group)
  map = node_table2[,2:3] %>% distinct(group, .keep_all = TRUE)
  mytable = node_table2 %>% df_mat(name,group,value) %>% as.matrix()
  mytable[is.na(mytable)] <- 0

  #--去除少于n个OTU的模块

  if (colSums(mytable)[colSums(mytable)< n] %>% length() == 0) {
    mytable_kp = mytable
  } else {
    mytable_kp = mytable[,-which(colSums(mytable)<n)]
  }


  head(mytable_kp)
  # 对应的处理一下map文件
  head(map)
  map_kp = map %>% filter(group %in% colnames(mytable_kp))
  mytable_kp = as.data.frame(mytable_kp)

  #---构造分组两两组合
  id = unique(node_table2$Group)

  network_pair = combn(id,2) %>% t() %>%
    as.matrix()


  total_mod_pairs = matrix(NA, nrow=nrow(network_pair), ncol=3)
  # i = 1
  for (i in 1:nrow(network_pair)){
    # 对两个组的每一个模块都进行比对
    module_pair = as.matrix(expand.grid(
      map_kp$group[which(map_kp$Group==network_pair[i,1])],
      map_kp$group[which(map_kp $Group==network_pair[i,2])]))

    total_mod_pairs[i,] = c(network_pair[i,], nrow(module_pair))
  }

  sig_mod_pairs = matrix(NA, nrow=0, ncol=4)
  sig_detailed_table = c("module1", "module2", "both", "P1A2", "P2A1", "A1A2", "p_raw", "p_adj")
  pairs_sum<-0
  i = 5
  for (i in 1:nrow(network_pair)){
    # 全部的需要比对的模块
    module_pair = as.matrix(expand.grid(map_kp$group[which(map_kp$Group==network_pair[i,1])],
                                        map_kp$group[which(map_kp $Group==network_pair[i,2])]))
    overlap = apply(module_pair, 1, FUN= find_overlap, bigtable= mytable_kp)
    only1 = apply(module_pair, 1, FUN= find_only_in_1, bigtable= mytable_kp)
    only2 = apply(module_pair, 1, FUN= find_only_in_2, bigtable= mytable_kp)
    denominator = apply(module_pair, 1, FUN= find_N, mapping=map, bigtable= mytable)
    none = denominator-(overlap + only1 + only2)
    count_table = data.frame(module1 = module_pair[,1],
                             module2 = module_pair[,2], Both=overlap, P1A2=only1, P2A1=only2, A1A2=none)

    pairs_num<-nrow(count_table)
    {
      pairs_sum<-sum(pairs_sum,pairs_num)
    }
    p_raw=c()
    # tt = 1
    for (tt in 1:nrow(count_table))
    {
      x=count_table[tt,]
      p = fisher_test(x)
      p_raw = c(p_raw, p)
    }

    count_table$p_raw = p_raw

    if (padj){
      count_table$p_adj = p.adjust(count_table$p_raw, method = "bonferroni")
    }else{
      count_table$p_adj = count_table$p_raw
    }

    network1 = network_pair[i,1]
    network2 = network_pair[i,2]
    sig_count = sum(count_table$p_adj<=0.05)
    # count_table$p_adj[3] = 0.001
    if(sig_count>0){
      sig_pairs_table = count_table[which(count_table$p_adj<=0.05),c(1:2)]
      sig_pairs_linked = paste(sig_pairs_table[,1], "-", sig_pairs_table[,2], sep="")
      sig_pairs = paste(sig_pairs_linked, collapse=",")

      sig_pairs_count_table = count_table[which(count_table$p_adj<=0.05),]
      row.names(sig_pairs_count_table) = sig_pairs_linked
    }else{
      sig_pairs = "None"
    }

    add_one_row = c(network1, network2, sig_count, sig_pairs)
    sig_mod_pairs = rbind(sig_mod_pairs, add_one_row)

    if(sig_count>0) sig_detailed_table = rbind(sig_detailed_table, sig_pairs_count_table)
    if(sig_count==0) sig_detailed_table = sig_detailed_table
    cat(i,": Module pairs from ",network_pair[i,]," have been calculated!!!\n")
  }
  sig_detailed_table<-sig_detailed_table[-1,]
  sig_detailed_table<-sig_detailed_table[!duplicated(sig_detailed_table),]
  pairs_num_pre<-nrow(sig_detailed_table)
  cat("\n","Totally ",pairs_sum," module pairs were tested.",sep = "")
  cat("\n","Totally ",pairs_num_pre," module pairs were preserved thanks to significant difference, accounted for ",
      paste(round(pairs_num_pre/pairs_sum*100,2),"%.",sep = ""),sep = "")

  file2=paste(export_path, "/Preserved module pairs (remove module(s) who owned less than ",n," node(s))",".txt",sep = "" )
  write.table(sig_detailed_table,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  return(list(data=sig_detailed_table,pairs_sum=pairs_sum,pairs_num_pre=pairs_num_pre))
}

"network.efficiency" <- function(graph){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  dd <- 1/shortest.paths(graph)
  diag(dd) <- NA
  efficiency <- mean(dd, na.rm=T)
  #denom <- nrow(dd)*(ncol(dd)-1)
  #sum(dd, na.rm=T)/denom
  return(efficiency)
}

"info.centrality.vertex" <- function(graph, net=NULL, verbose=F){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  if(is.null(net)) net <- network.efficiency(graph)
  if(is.numeric(net)==F){
    warning("Please ensure net is a scalar numeric")
    net <- network.efficiency(graph)
  }
  count <- c()
  for(i in 1:length(V(graph))){
    count <- c(count, (net-network.efficiency(delete.vertices(graph, i)))/net)
    if(verbose){
      print(paste("node",i,"current\ info\ score", count[i], collapse="\t"))
    }
  }
  return(count)
}

"info.centrality.network" <- function(graph, net=network.efficiency(graph), verbose=F) sum(info.centrality.vertex(graph))



"unique.char" <- function(sas) {
  all.let<-strsplit(sas, NULL)
  chac<-c()
  for (tt in 1:length(all.let)) {
    cha<-all.let[[tt]]
    {
      chac<-c(chac,cha)
    }
  }
  chac<-unique(chac)
  return(chac)

}

"kruskalmuticomp.t" <- function(input.data) {
  input.data1<-data.frame()
  for (index.n in unique(input.data$variable)) {
    input.data.new<-input.data[which(input.data$variable==index.n),]
    input.datax.new<-kruskalmuticomp(input.data.new,'variable','group','value',
                                     reversed = FALSE)

    input.data1.new<-input.datax.new
    match.n<-c(letters[1:5],"ab","bc","cd","de",
               "abc","bcd","cde",
               "abcd","bcde",
               "abced")
    names(match.n)<-c(10000,20000,30000,40000,50000,
                      12000,23000,34000,45000,
                      12300,23400,34500,
                      12340,23450,
                      12345)
    match.n<-match.n%>%data.frame()
    match.n$value<-rownames(match.n)
    match.n$value<-match.n$value%>%as.numeric()
    colnames(match.n)[1]<-"letter"
    input.data1.new$Letter<-match.n$value[match(input.data1.new$Letters,match.n$letter)]
    max.mean<-which(input.data1.new$mean==max(input.data1.new$mean))[1]
    max.let<-which(input.data1.new$Letter==max(input.data1.new$Letter))[1]
    if (max.mean==max.let) input.datax.new<-input.datax.new
    if (max.mean!=max.let) input.datax.new<-kruskalmuticomp(input.data.new,'variable','group','value',
                                                            reversed = TRUE)

    {
      input.data1<-rbind(input.data1,input.datax.new)
    }
  }
  return(input.data1)
}

"keep.decmi" <- function(x,decmi=3) {
  tt<-sprintf(paste("%0.",decmi,"f",sep = ""), round(x,decmi))
  return(tt)
}

"remove.letter" <- function(x)
{
  sapply(
    lapply(strsplit(x, NULL), function(y){
      if (length(which(y %in% letters)) !=0) y<-y[-which(y %in% letters)] else y<-y
      if (length(which(y %in% " ")) !=0) y<-y[-which(y %in% " ")] else y<-y
    }),
    paste, collapse = "")
}

"keep.letter" <- function(x)
{
  sapply(
    lapply(strsplit(x, NULL), function(y){
      if (length(which(y %in% letters)) !=0) y<-y[which(y %in% letters)] else y<-""
    }),
    paste, collapse = "")
}

"install.allpackage" <- function() {
  x<-c("vegan", "dplyr", "BiocManager","ggrepel", "doParallel", "picante","vegan","ggpubr",
       "devtools","multcomp","multcompView","pgirmess","rstatix","rdacca.hp",
       "patchwork","ggVennDiagram","ape","ggrepel","circlize","tidyverse",
       "cowplot","dplyr","doParallel","reshape2","UpSetR","doBy","Hmisc",
       "ggtext","plotrix","stringr","foreach" ,"mgcv", "reshape2", "ggplot2")

  x <- x[!x %in% installed.packages()]
  lapply(x, install.packages, character.only = TRUE)

  ### devtools
  y<-c("xiangpin/ggtreeExtra",
       "YuLab-SMU/ggtree",
       "YuLab-SMU/MicrobiotaProcess",
       "YuLab-SMU/treeio",
       "GuillemSalazar/EcolUtils",
       "pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  z<-c("phyloseq",
       "edgeR",
       "limma")

  aa<-strsplit(y,"/")
  y.n<-c()
  for (tt in 1:length(aa)) {
    sa<-aa[[tt]][length(aa[[tt]])]
    {y.n<-c(y.n,sa)}
  }
  y.nn <- y.n[!y.n %in% installed.packages()]
  y<-y[y.nn==y.n]
  z <- z[!z %in% installed.packages()]
  require(devtools)
  require(BiocManager)
  lapply(y, devtools::install_github, force = TRUE)
  BiocManager::install(z,force = TRUE)
}


"calc_indivGroup_mean" <- function(otutab) {
  colname<-colnames(otutab)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)
  group_num<-length(trt_id)
  groupname<-substr(colnames(otutab),start = 1,stop = 2)%>%unique()
  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otutab))
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
          grep(x,colnames(otutab))
        }),2,FUN = function(x){
          otutab[,x]
        }),
      function(x){
        rowMeans(x)
      })
  } else {
    split_otu <- lapply(
      lapply(
        sapply(trt_id,function(x){
          grep(x,colnames(otutab))
        }),FUN = function(x){
          otutab[,x]
        }),
      function(x){
        rowMeans(x)
      })
  }
  return(split_otu)

}

"calc_indivGroup_abun" <- function(otutab) {
  colname<-colnames(otutab)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)
  group_num<-length(trt_id)
  groupname<-substr(colnames(otutab),start = 1,stop = 2)%>%unique()
  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otutab))
  }),function(y){
    y %>%length()
  })
  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()

  if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otux <- apply(
      sapply(trt_id,function(x){
        grep(x,colnames(otutab))
      }),2,FUN = function(x){
        otutab[,x]
      })
  } else {
    split_otux <- lapply(
      sapply(trt_id,function(x){
        grep(x,colnames(otutab))
      }),FUN = function(x){
        otutab[,x]
      })
  }
  return(split_otux)
}
