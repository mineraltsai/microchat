"circosfile" <- function(df) {
  for (t in 1:length(colnames(df))) {
    firstname<-tolower(substr(colnames(df[t]),1,1))
    df[,t][which(df[,t]=="")]<-paste(firstname,"__unclassified",sep="")
  }

  for (k in 2:length(colnames(df))) {
    firstname<-tolower(substr(colnames(df[k]),1,1))
    if (length(which(df[,k]==paste(firstname,"__norank",sep=""))) !=0) {
      duppos<-which(df[,k]==paste(firstname,"__norank",sep=""))

      if (!is.null(duppos)) {
        tto<-df[duppos,k-1]
        newtto<-cbind(tto,duppos)%>%data.frame()
        newtto$duppos<-newtto$duppos%>%as.numeric()
        newtto<-newtto[order(newtto$tto),]
        dupnum<-length(unique(newtto$tto))

        if (dupnum ==1) {
          df[duppos,k]=paste(firstname,"__norank",sep="")
        } else {

          for (tt in 1:dupnum) {
            kpos<-newtto[which(newtto$tto == unique(newtto$tto)[tt]),]$duppos
            df[kpos,k]<-paste(df[kpos,k],tt,sep = "")
          }


        }


      }

    }


  }

  for (k in 2:length(colnames(df))) {
    firstname<-tolower(substr(colnames(df[k]),1,1))
    if (length(which(df[,k]==paste(firstname,"__unclassified",sep=""))) !=0) {
      duppos<-which(df[,k]==paste(firstname,"__unclassified",sep=""))

      if (!is.null(duppos)) {
        tto<-df[duppos,k-1]
        newtto<-cbind(tto,duppos)%>%data.frame()
        newtto$duppos<-newtto$duppos%>%as.numeric()
        newtto<-newtto[order(newtto$tto),]
        dupnum<-length(unique(newtto$tto))

        if (dupnum ==1) {
          df[duppos,k]=paste(firstname,"__unclassified",sep="")
        } else {

          for (tt in 1:dupnum) {
            kpos<-newtto[which(newtto$tto == unique(newtto$tto)[tt]),]$duppos
            df[kpos,k]<-paste(df[kpos,k],tt,sep = "")
          }

        }

      }

    }

  }

  for (k in 2:length(colnames(df))) {
    firstname<-tolower(substr(colnames(df[k]),1,1))
    if (length(which(substr(df[,k],start = 1,stop = 5) ==paste(firstname,"__un",sep = "")))!=0) {
      duppos<-which(substr(df[,k],start = 1,stop = 5) ==paste(firstname,"__un",sep = ""))

      dd<-which(substr(df[,k],start = 1,stop = 10) ==paste(firstname,"__unclass",sep = ""))

      duppos1<-setdiff(duppos,dd)
      duppos<-duppos1
      if (!is.null(duppos)) {
        tto<-df[duppos,k-1]
        newtto<-cbind(tto,duppos)%>%data.frame()
        newtto$duppos<-newtto$duppos%>%as.numeric()
        newtto<-newtto[order(newtto$tto),]

        dupnum<-length(unique(newtto$tto))

        if (dupnum !=1) {

          for (tt in 1:dupnum) {
            kpos<-newtto[which(newtto$tto == unique(newtto$tto)[tt]),]$duppos
            df[kpos,k]<-paste(df[kpos,k],tt,sep = "")
          }

        }

      }

    }

  }
  return(df)
}

"circosfile_compos" <- function(ddf,taxon,edge_circos,export_path) {
  options(warn = -1)
  allcolor<-taxon$Phylum%>%unique()%>%data.frame()
  allcolor$color<-paste("silver",1:nrow(allcolor),sep = "")
  colnames(allcolor)<-c("taxa","color")

  Phylumnum<-ddf[,3:4]
  Phylumfreq<-data.frame(Phylumnum$Phylum %>% table())
  colnames(Phylumfreq)<-c('Phylum','Freq')

  Phylumsort <- arrange(Phylumfreq,Phylumfreq$Phylum,.by_group = TRUE)
  cml<-merge(Phylumsort,Phylumfreq,by = 'Phylum')
  cumsum(cml$Freq.x)

  colnames(cml)<- c('Phylum','start',"end")

  hsPhylum<-cml
  hsPhylum$color1<-"fill_color="

  hsPhylum$color2<-"silver"
  hsPhylum$nodenum<-c(1:length(hsPhylum$Phylum))
  silver<-hsPhylum
  hsPhylum<-unite(hsPhylum,'rgb',color1,color2,nodenum,sep = "")
  hsPhylum$plus<-cumsum(hsPhylum$end)
  hsPhylum
  hsPhylum$end<-hsPhylum$plus*2
  hsPhylum$plus <- hsPhylum$end
  hsPhylum$start[1]<-0
  hsPhylum$start[2:length(hsPhylum$Phylum)]<-hsPhylum$plus[1:(length(hsPhylum$Phylum)-1)]
  classcircle<-hsPhylum
  hsPhylum$chr<-"k__Acidobacteria"
  hstrackouter<-hsPhylum

  hstrackouter<-dplyr::select(hstrackouter,6,2,3,1,dplyr::everything())

  hstrackouter<-hstrackouter[,-5]
  hstrackouter<-hstrackouter[,-5]
  hstrackouter<-unite(hstrackouter,"outertrack", sep = " ")

  flie=paste(export_path, "outertrack",".txt", sep = "")
  write.table(hstrackouter,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  hsPhylum$band<-"band"
  hsPhylum$Phylum1<-hsPhylum$Phylum
  hsPhylum$color2<-"silver"
  hsPhylum$nodenum<-c(1:length(hsPhylum$Phylum))
  hsPhylum<-unite(hsPhylum,'rgb',color2,nodenum,sep = "")
  hsPhylum<-left_join(hsPhylum,allcolor,by=c("Phylum"="taxa"))

  hsPhylum1<-dplyr::select(hsPhylum,6,5,1,7,2,3,9,dplyr::everything())
  hsPhylum1<-hsPhylum1[,-(8:9)]
  maxkingdom<-max(hsPhylum1$end)
  tt<-c("chr", '-', "k__Acidobacteria", "k__Acidobacteria", 0, maxkingdom, "silver1")
  ttt<-t(tt)%>%data.frame()
  colnames(ttt)<-colnames(hsPhylum1)

  hsPhylum3<-rbind(ttt,hsPhylum1)

  hsPhylum2<-unite(hsPhylum1,'colors',sep = " ")
  hskingdom2<-unite(hsPhylum3,'colors',sep = " ")

  flie=paste(export_path, "hsphylum2",".txt", sep = "")
  write.table(hsPhylum2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  flie=paste(export_path, "kingdom",".txt", sep = "")
  write.table(hskingdom2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  ##########-------------


  Classnum<-ddf[,3:4]


  phycla<-unite(Classnum,"phycla",Phylum,Class, sep = "" )
  phycla1<-data.frame(phycla %>% table())
  colnames(phycla1)<-c('Class','Freq')


  phycla2<-filter(phycla1,phycla1$Freq>0)

  Classsort <- arrange(phycla2,phycla2$Class,.by_group = TRUE)
  cml<-Classsort
  cumsum(cml$Freq)
  colnames(cml)<- c('Class','end')

  hsClass<-cml
  hsClass$color1<-"fill_color="

  hsClass$color2<-"silver"
  hsClass$nodenum<-c(1:length(hsClass$Class))
  silver<-hsClass
  hsClass<-unite(hsClass,'rgb',color1,color2,nodenum,sep = "")
  hsClass$plus<-cumsum(hsClass$end)
  hsClass
  hsClass$end<-hsClass$plus*2
  hsClass$plus <- hsClass$end
  hsClass$start[1]<-0
  hsClass$start[2:length(hsClass$Class)]<-hsClass$plus[1:(length(hsClass$Class)-1)]
  classcircle<-hsClass

  {
    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$Class[t],start = 1,stop = 10) ,
                      substr(hsPhylum1$Phylum,start = 1,stop = 10))
      classcircle$color[t]<-hsPhylum1$color[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color


    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")

  }
  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))
  classcircle2<-unite(classcircle2,'colors',sep = " ")
  flie=paste(export_path, "hsclass2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  ########------
  ordernum<-ddf[,3:5]

  phycla<-unite(ordernum,"phycla",Phylum,Class,Order, sep = "" )
  phycla1<-data.frame(phycla %>% table())
  colnames(phycla1)<-c('Order','Freq')

  phycla2<-filter(phycla1,phycla1$Freq>0)

  ordersort <- arrange(phycla2,phycla2$Order,.by_group = TRUE)
  cml<-ordersort
  cumsum(cml$Freq)

  colnames(cml)<- c('Order',"end")

  hsorder<-cml
  hsorder$color1<-"fill_color="

  hsorder$color2<-"silver"
  hsorder$nodenum<-c(1:length(hsorder$Order))
  silver<-hsorder
  hsorder<-unite(hsorder,'rgb',color1,color2,nodenum,sep = "")
  hsorder$plus<-cumsum(hsorder$end)
  hsorder
  hsorder$end<-hsorder$plus*2
  hsorder$plus <- hsorder$end
  hsorder$start[1]<-0
  hsorder$start[2:length(hsorder$Order)]<-hsorder$plus[1:(length(hsorder$Order)-1)]
  classcircle<-hsorder

  {
    classcircle$Order<-classcircle$Order%>%as.character()
    rrr<-strsplit(classcircle$Order, "o__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1
    stopnum<-str_length(classcircle$pretax)%>%max()
    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Class,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")



  }
  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))
  classcircle2<-unite(classcircle2,'colors',sep = " ")
  flie=paste(export_path, "hsorder2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)



  ########------
  Familynum<-ddf[,3:6]


  phycla<-unite(Familynum,"phycla",Phylum,Class,Order,Family, sep = "" )
  phycla1<-data.frame(phycla %>% table())
  colnames(phycla1)<-c('Family','Freq')

  phycla2<-filter(phycla1,phycla1$Freq>0)
  Familysort <- arrange(phycla2,phycla2$Family,.by_group = TRUE)
  cml<-Familysort
  cumsum(cml$Freq)

  colnames(cml)<- c('Family',"end")

  hsFamily<-cml
  hsFamily$color1<-"fill_color="

  hsFamily$color2<-"silver"
  hsFamily$nodenum<-c(1:length(hsFamily$Family))
  silver<-hsFamily
  hsFamily<-unite(hsFamily,'rgb',color1,color2,nodenum,sep = "")
  hsFamily$plus<-cumsum(hsFamily$end)
  hsFamily
  hsFamily$end<-hsFamily$plus*2
  hsFamily$plus <- hsFamily$end
  hsFamily$start[1]<-0
  hsFamily$start[2:length(hsFamily$Family)]<-hsFamily$plus[1:(length(hsFamily$Family)-1)]
  classcircle<-hsFamily

  {
    classcircle$Family<-classcircle$Family%>%as.character()
    rrr<-strsplit(classcircle$Family, "f__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1

    stopnum<-str_length(classcircle$pretax)%>%max()

    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Order,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")


  }


  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))

  classcircle2<-unite(classcircle2,'colors',sep = " ")

  flie=paste(export_path, "hsfamily2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  ########------
  Genusnum<-ddf[,3:7]


  phycla<-unite(Genusnum,"phycla",Phylum,Class,Order,Family,Genus, sep = "" )
  phycla1<-data.frame(phycla %>% table())
  colnames(phycla1)<-c('Genus','Freq')

  phycla2<-filter(phycla1,phycla1$Freq>0)
  Genussort <- arrange(phycla2,phycla2$Genus,.by_group = TRUE)
  cml<-Genussort
  cumsum(cml$Freq)

  colnames(cml)<- c('Genus',"end")

  hsGenus<-cml
  hsGenus$color1<-"fill_color="

  hsGenus$color2<-"silver"
  hsGenus$nodenum<-c(1:length(hsGenus$Genus))
  silver<-hsGenus
  hsGenus<-unite(hsGenus,'rgb',color1,color2,nodenum,sep = "")
  hsGenus$plus<-cumsum(hsGenus$end)
  hsGenus
  hsGenus$end<-hsGenus$plus*2
  hsGenus$plus <- hsGenus$end
  hsGenus$start[1]<-0
  hsGenus$start[2:length(hsGenus$Genus)]<-hsGenus$plus[1:(length(hsGenus$Genus)-1)]
  classcircle<-hsGenus

  {
    classcircle$Genus<-classcircle$Genus%>%as.character()
    rrr<-strsplit(classcircle$Genus, "g__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1

    stopnum<-str_length(classcircle$pretax)%>%max()

    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Family,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")


  }


  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))
  classcircle2<-unite(classcircle2,'colors',sep = " ")

  flie=paste(export_path, "hsgenus2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  ########------
  Speciesnum<-subset(ddf,select = c(3:8,1))
  genusgene<-unite(ddf,"phy",Phylum,Class,Order,Family,Genus,Species,geneid, sep = "" )


  phycla<-unite(Speciesnum,"phycla",Phylum,Class,Order,Family,Genus,Species,geneid, sep = "" )
  phycla1<-data.frame(phycla %>% table())
  colnames(phycla1)<-c('Species','Freq')

  phycla2<-filter(phycla1,phycla1$Freq>0)

  Genussort <- arrange(phycla2,phycla2$Species,.by_group = TRUE)
  cml<-Genussort
  cumsum(cml$Freq)

  colnames(cml)<- c('Species',"end")

  hsGenus<-cml
  hsGenus$color1<-"fill_color="

  hsGenus$color2<-"silver"
  hsGenus$nodenum<-c(1:length(hsGenus$Species))
  silver<-hsGenus
  hsGenus<-unite(hsGenus,'rgb',color1,color2,nodenum,sep = "")
  hsGenus$plus<-cumsum(hsGenus$end)
  hsGenus
  hsGenus$end<-hsGenus$plus*2
  hsGenus$plus <- hsGenus$end
  hsGenus$start[1]<-0
  hsGenus$start[2:length(hsGenus$Species)]<-hsGenus$plus[1:(length(hsGenus$Species)-1)]
  classcircle<-hsGenus

  {
    classcircle$Species<-classcircle$Species%>%as.character()
    rrr<-strsplit(classcircle$Species, "s__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1

    stopnum<-str_length(classcircle$pretax)%>%max()

    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Genus,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")


  }

  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  hsGenus<-classcircle1
  hsGenus$phy<-hsGenus$Species
  genus<-merge(hsGenus,genusgene,by="phy")
  hsspe<-genus
  hsspe<-cbind(hsspe,Speciesnum$geneid)

  hsspe<-subset(hsspe,select=c(11,6,3,15))
  hsspe1<-unite(hsspe,'colors',sep = " ")
  flie=paste(export_path, "hsotu2",".txt", sep = "")
  write.table(hsspe1,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  ########------

  colnames(hsspe)<-c('chr','start','end',"geneid")
  scatterotu<-merge(ddf,hsspe,by="geneid")
  scatterotu<-scatterotu[,9:12]
  scatterotu1<-subset(scatterotu,select=c(2,3,4,1))
  scatterotu1$fc<-scatterotu1$fc%>%as.numeric()
  scatterotu1$fc<-log2(scatterotu1$fc)
  scatterotu2<-unite(scatterotu1,'colors',sep = " ")
  flie=paste(export_path, "scatterotu2",".txt", sep = "")
  write.table(scatterotu2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  hsgeneid<-subset(genus,select=c(chr,start,end,newcolor1))
  hsgeneid2<-unite(hsgeneid,'colors',sep = " ")
  flie=paste(export_path, "hsgeneid2",".txt", sep = "")
  write.table(hsgeneid2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  hsspe
  edge = edge_circos
  edge1<-edge
  hsotuedge<-hsspe
  hsotuedge$node2<-hsotuedge$geneid
  colnames(hsotuedge)<-c('chr','start','end','node1','node2')
  colnames(edge1)<-c('node1',"node2",'cor')
  edgestart<-merge(edge1,hsotuedge, by ="node1")
  edgestart<-edgestart[,-7]
  colnames(edgestart)<-c('node-start1','node2','cor','chr1','start1','start2')
  edgetotal<-edgestart
  edges<-merge(edgetotal,hsotuedge, by ="node2")
  edges<-edges[,-1]
  edges<-edges[,-1]
  edges<-edges[,-8]
  edges<-subset(edges,select=c(2,3,4,5,6,7,1))

  edge1<-filter(edges,edges$cor>0)
  edge1<-edge1[,-7]
  edge1<-unite(edge1,'edge',sep = " ")

  flie=paste(export_path, "linkpos",".txt", sep = "")
  write.table(edge1,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  edge2<-filter(edges,edges$cor<0)
  edge2<-edge2[,-7]
  edge2<-unite(edge2,'edge',sep = " ")

  flie=paste(export_path, "linkneg",".txt", sep = "")
  write.table(edge2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

}

"circosfile_abun" <- function(ddf,taxon,edge_circos,export_path) {
  options(warn = -1)
  ddf$fc<-ddf$fc%>%as.numeric()

  allcolor<-taxon$Phylum%>%unique()%>%data.frame()
  allcolor$color<-paste("silver",1:nrow(allcolor),sep = "")
  colnames(allcolor)<-c("taxa","color")



  Phylumnum<-subset(ddf,select=c(3,9))

  ###unite

  Phylumnum<-Phylumnum%>%
    group_by(Phylum)%>%
    summarise_all(sum)

  Phylumnum$start<-Phylumnum$fc
  Phylumnum$end<-cumsum(Phylumnum$fc)
  Phylumnum$end<-Phylumnum$end*2
  Phylumnum$start[1]<-0
  Phylumnum$start[2:length(Phylumnum$Phylum)]<-Phylumnum$end[1:(length(Phylumnum$Phylum)-1)]

  hsPhylum<-Phylumnum
  classcircle<-hsPhylum


  hsPhylum$chr<-"k__Acidobacteria"
  hstrackouter<-hsPhylum

  hstrackouter<-subset(hstrackouter,select=c(5,3,4,1))
  hstrackouter<-unite(hstrackouter,"outertrack", sep = " ")

  flie=paste(export_path, "outertrack",".txt", sep = "")
  write.table(hstrackouter,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  hsPhylum$band<-"band"
  hsPhylum$Phylum1<-hsPhylum$Phylum
  hsPhylum$color2<-"silver"
  hsPhylum$nodenum<-c(1:length(hsPhylum$Phylum))
  hsPhylum<-unite(hsPhylum,'rgb',color2,nodenum,sep = "")

  hsPhylum<-left_join(hsPhylum,allcolor,by=c("Phylum"="taxa"))

  hsPhylum1<-subset(hsPhylum,select=c(6,5,1,7,3,4,9))
  maxkingdom<-max(hsPhylum1$end)
  tt<-c("chr", '-', "k__Acidobacteria", "k__Acidobacteria", 0, maxkingdom, "silver1")
  ttt<-t(tt)%>%data.frame()
  colnames(ttt)<-colnames(hsPhylum1)

  hsPhylum3<-rbind(ttt,hsPhylum1)

  hsPhylum2<-unite(hsPhylum1,'colors',sep = " ")
  hskingdom2<-unite(hsPhylum3,'colors',sep = " ")

  flie=paste(export_path, "hsphylum2",".txt", sep = "")
  write.table(hsPhylum2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  flie=paste(export_path, "kingdom",".txt", sep = "")
  write.table(hskingdom2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  ##########-------------
  Classnum<-subset(ddf,select=c(3:4,9))

  classcla<-unite(Classnum,"Phylum",Phylum,Class, sep = "" )

  classcla<-classcla%>%
    group_by(Phylum)%>%
    summarise_all(sum)

  classcla$start<-classcla$fc
  classcla$end<-cumsum(classcla$fc)
  classcla$end<-classcla$end*2
  classcla$start[1]<-0
  classcla$start[2:length(classcla$Phylum)]<-classcla$end[1:(length(classcla$Phylum)-1)]

  hsPhylum<-classcla
  classcircle<-classcla
  classcircle$color<-1
  {

    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$Phylum[t],start = 1,stop = 10) ,
                      substr(hsPhylum1$Phylum,start = 1,stop = 10))
      classcircle$color[t]<-hsPhylum1$color[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color


    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")

  }
  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))

  classcircle2<-unite(classcircle2,'colors',sep = " ")
  flie=paste(export_path, "hsclass2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)



  ########------
  ordernum<-subset(ddf,select=c(3:5,9))

  classcla<-unite(ordernum,"Phylum",Phylum,Class,Order, sep = "" )

  classcla<-classcla%>%
    group_by(Phylum)%>%
    summarise_all(sum)

  classcla$start<-classcla$fc
  classcla$end<-cumsum(classcla$fc)
  classcla$end<-classcla$end*2
  classcla$start[1]<-0
  classcla$start[2:length(classcla$Phylum)]<-classcla$end[1:(length(classcla$Phylum)-1)]

  hsPhylum<-classcla
  classcircle<-classcla
  classcircle$color<-1

  {
    classcircle$Phylum<-classcircle$Phylum%>%as.character()
    rrr<-strsplit(classcircle$Phylum, "o__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1
    stopnum<-str_length(classcircle$pretax)%>%max()
    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Phylum,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")



  }
  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))

  classcircle2<-unite(classcircle2,'colors',sep = " ")
  flie=paste(export_path, "hsorder2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  ########------
  familynum<-subset(ddf,select=c(3:6,9))

  classcla<-unite(familynum,"Phylum",Phylum,Class,Order,Family, sep = "" )

  classcla<-classcla%>%
    group_by(Phylum)%>%
    summarise_all(sum)

  classcla$start<-classcla$fc
  classcla$end<-cumsum(classcla$fc)
  classcla$end<-classcla$end*2
  classcla$start[1]<-0
  classcla$start[2:length(classcla$Phylum)]<-classcla$end[1:(length(classcla$Phylum)-1)]

  hsPhylum<-classcla
  classcircle<-classcla
  classcircle$color<-1


  {
    classcircle$Phylum<-classcircle$Phylum%>%as.character()
    rrr<-strsplit(classcircle$Phylum, "f__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1

    stopnum<-str_length(classcircle$pretax)%>%max()

    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Phylum,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")


  }


  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))

  classcircle2<-unite(classcircle2,'colors',sep = " ")

  flie=paste(export_path, "hsfamily2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)



  ########------
  genusnum<-subset(ddf,select=c(3:7,9))

  classcla<-unite(genusnum,"Phylum",Phylum,Class,Order,Family,Genus, sep = "" )

  classcla<-classcla%>%
    group_by(Phylum)%>%
    summarise_all(sum)

  classcla$start<-classcla$fc
  classcla$end<-cumsum(classcla$fc)
  classcla$end<-classcla$end*2
  classcla$start[1]<-0
  classcla$start[2:length(classcla$Phylum)]<-classcla$end[1:(length(classcla$Phylum)-1)]

  hsPhylum<-classcla
  classcircle<-classcla
  classcircle$color<-1

  {
    classcircle$Phylum<-classcircle$Phylum%>%as.character()
    rrr<-strsplit(classcircle$Phylum, "g__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1

    stopnum<-str_length(classcircle$pretax)%>%max()

    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Phylum,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")


  }


  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  classcircle2<-subset(classcircle1,select=c(chr,start,end,newcolor1))

  classcircle2<-unite(classcircle2,'colors',sep = " ")

  flie=paste(export_path, "hsgenus2",".txt", sep = "")
  write.table(classcircle2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)



  ########------

  specnum<-subset(ddf,select=c(3:8,1,9))

  classcla<-unite(specnum,"Phylum",Phylum,Class,Order,Family,Genus,Species,geneid, sep = "" )

  classcla<-classcla%>%
    group_by(Phylum)%>%
    summarise_all(sum)

  classcla$start<-classcla$fc
  classcla$end<-cumsum(classcla$fc)
  classcla$end<-classcla$end*2
  classcla$start[1]<-0
  classcla$start[2:length(classcla$Phylum)]<-classcla$end[1:(length(classcla$Phylum)-1)]

  hsPhylum<-classcla
  classcircle<-classcla
  classcircle$color<-1

  {
    classcircle$Phylum<-classcircle$Phylum%>%as.character()
    rrr<-strsplit(classcircle$Phylum, "s__",fixed = TRUE)%>%data.frame()%>%t()%>%data.frame()
    classcircle$pretax<-rrr$X1

    stopnum<-str_length(classcircle$pretax)%>%max()

    for (t in 1:nrow(classcircle)) {
      matchnum<-match(substr(classcircle$pretax[t],start = 1,stop = stopnum) ,
                      substr(classcircle1$Phylum,start = 1,stop = stopnum))
      classcircle$color[t]<-classcircle1$newcolor[matchnum]

    }

    classcircle$error<-classcircle$end-classcircle$start
    classcircle$newcolor<-classcircle$color

    for (kk in 1:length(unique(classcircle$color))) {
      hhh<-which(classcircle$color==unique(classcircle$color)[kk])
      if (length(hhh)==1) {
        classcircle$newcolor[hhh]<-classcircle$color[hhh]
      } else {
        relpos<-which(classcircle[hhh,]$error==max(classcircle[hhh,]$error))[1]
        abspos<-hhh[relpos]
        otherpos<-setdiff(hhh,abspos)
        classcircle$newcolor[abspos]<-classcircle$color[abspos]
        classcircle$newcolor[otherpos]<-""
      }
    }

    maxnumcal<-classcircle$newcolor[which(classcircle$newcolor!="")]
    maxnum<-substr(maxnumcal,start = 7,stop = 10)%>%as.numeric()%>%max()
    nullcolnum<-classcircle[which(classcircle$newcolor==""),]%>%nrow()
    classcircle[which(classcircle$newcolor==""),]$newcolor<-paste("silver",(maxnum+1):(maxnum+nullcolnum),sep = "")


  }

  classcircle1<-classcircle
  classcircle1$chr<-"k__Acidobacteria"
  classcircle1$newcolor1<-paste(rep("fill_color=",nrow(classcircle1)),
                                classcircle1$newcolor,sep = "")

  hsGenus<-classcircle1
  hsspe<-hsGenus
  hsspe<-cbind(hsspe,specnum$geneid)

  hsspe<-subset(hsspe,select=c(chr,start,end,11))
  hsspe1<-unite(hsspe,'colors',sep = " ")
  flie=paste(export_path, "hsotu2",".txt", sep = "")
  write.table(hsspe1,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  ########------

  colnames(hsspe)<-c('chr','start','end',"geneid")
  scatterotu<-merge(ddf,hsspe,by="geneid")
  scatterotu<-scatterotu[,9:12]
  scatterotu1<-subset(scatterotu,select=c(2,3,4,1))
  scatterotu1$fc<-scatterotu1$fc%>%as.numeric()
  scatterotu1$fc<-log2(scatterotu1$fc)
  scatterotu2<-unite(scatterotu1,'colors',sep = " ")
  flie=paste(export_path, "scatterotu2",".txt", sep = "")
  write.table(scatterotu2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  hsgeneid<-subset(hsGenus,select=c(chr,start,end,newcolor1))

  hsgeneid2<-unite(hsgeneid,'colors',sep = " ")
  flie=paste(export_path, "hsgeneid2",".txt", sep = "")
  write.table(hsgeneid2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


  hsspe
  edge = edge_circos
  edge1<-edge
  hsotuedge<-hsspe
  hsotuedge$node2<-hsotuedge$geneid
  colnames(hsotuedge)<-c('chr','start','end','node1','node2')
  colnames(edge1)<-c('node1',"node2",'cor')
  edgestart<-merge(edge1,hsotuedge, by ="node1")
  edgestart<-edgestart[,-7]
  colnames(edgestart)<-c('node-start1','node2','cor','chr1','start1','start2')
  edgetotal<-edgestart
  edges<-merge(edgetotal,hsotuedge, by ="node2")
  edges<-edges[,-1]
  edges<-edges[,-1]
  edges<-edges[,-8]
  edges<-subset(edges,select=c(2,3,4,5,6,7,1))


  edge1<-filter(edges,edges$cor>0)
  edge1<-edge1[,-7]
  edge1<-unite(edge1,'edge',sep = " ")

  flie=paste(export_path, "linkpos",".txt", sep = "")
  write.table(edge1,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)



  edge2<-filter(edges,edges$cor<0)
  edge2<-edge2[,-7]
  edge2<-unite(edge2,'edge',sep = " ")

  flie=paste(export_path, "linkneg",".txt", sep = "")
  write.table(edge2,file = flie, row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)


}

"rmnet_thres_calc" <- function(abun,
                               filter_num = 1,
                               seed=152510,
                               cor.method = "pearson") {
  library(igraph)
  library(dplyr)
  library(ggtext)
  library(tidyr)
  library(ggpubr)
  options(warn = -1)
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


  mat <- function(g) {
    gmat<-lapply(g, function(x){
      gg<-as.matrix(x)
      return(gg)
    })
    return(gmat)
  }
  ###subgroup
  otus<-mat(otu)

  nor <- lapply(otus,function(x){
    new<-scale(t(x))##
    occor.r = cor(new,method = cor.method)
    diag(occor.r) <- 0
    occor.r[is.na(occor.r)]=0
    return(occor.r)
  })

  trvalue <- lapply(nor,function(x){
    #Get threshold by RMT
    res <- rm.get.threshold(x,discard.zeros = TRUE,plot.spacing=FALSE,
                            save.fit = FALSE,
                            unfold.method = "gaussian")

    xas<-res$tested.thresholds
    yas<-res$p.ks
    plot<-cbind(xas,yas)%>%data.frame()
    threshold<-plot$xas[which(plot$yas == max(plot$yas))]
    thre <- threshold
    return(list(bor=x,plot=plot,result=thre))
  })

  return(list(trvalue=trvalue,split_otu=split_otu))
}

"rmnet_thres_select" <- function(trvalue,
                                 selected.thres=c(0.65,0.75,0.5,0.6,0.7),
                                 export_path="microbial network analysis") {
  export_path<-paste(export_path,"/microbial network analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  if (length(selected.thres) <=1) {

    st.trvalue<-lapply(trvalue, function(x){

      x$result<-selected.thres
      return(x)
    })

  } else {

    if (length(selected.thres)<length(trvalue)) {
      message("\n","selected.thres is not enough for the required number","\n")
    } else {

      for (tt in 1:length(trvalue)) {
        trvalue[[tt]]$result<-selected.thres[[tt]]
      }
      st.trvalue<-trvalue
    }
  }

  net <- lapply(st.trvalue,function(x){
    cleaned.matrix <- rm.denoise.mat(x$bor, threshold = x$result) #denoise adjacency matrix by threshold
    cleaned.matrix <- rm.discard.zeros(cleaned.matrix) #delete all-zero rows
    g = graph_from_adjacency_matrix(cleaned.matrix ,mode="undirected",weighted=TRUE,diag=FALSE)
    result<-x$result
    name<-names(x)
    return(list(g=g,cleaned.matrix=cleaned.matrix,result=result,name=name))
  })


  g<-lapply(net, function(x){
    gg<-x$g
    return(gg)
  })

  netmod<-data.frame()
  for(k in 1:length(g)){
    name<-names(g)[k]
    g1 <- g[[k]]

    E(g1)$correlation <- E(g1)$weight

    edgelen<-E(g1)$correlation
    alllen<-edgelen%>%length()
    neglen<-edgelen[edgelen<0]%>%length()
    poslen<-edgelen[edgelen>0]%>%length()

    pos.pro<-paste(round(poslen/alllen*100,2),"%")
    neg.pro<-paste(round(neglen/alllen*100,2),"%")

    pos<-paste(poslen," (",pos.pro,")",sep = "")
    neg<-paste(neglen," (",neg.pro,")",sep = "")

    E(g1)$weight <- abs(E(g1)$weight)
    set.seed(007)

    V(g1)$modularity <- membership(cluster_fast_greedy(g1))
    modulenum<-length(cluster_fast_greedy(g1))
    V(g1)$label <- V(g1)$name
    V(g1)$label <- NA
    modu_sort <- V(g1)$modularity %>% table() %>% sort(decreasing = T)
    node_table = as.data.frame(vertex.attributes(g1))
    node_table$degree = igraph::degree(g1)
    node_table$closeness=closeness(g1, mode="all")

    module_membership = membership(cluster_fast_greedy(g1))
    node_table$module= membership(cluster_fast_greedy(g1))
    edge_table = as.data.frame(edge.attributes(g1))

    netmod1<-cbind(name,modulenum,nrow(node_table),nrow(edge_table),pos,neg)
    colnames(netmod1)<-c("group","module number","node number","edge number","positive edge","negative edge")
    {
      netmod<-rbind(netmod,netmod1)
    }
  }

  result<-lapply(st.trvalue, function(x){
    result<-x$result
    return(result)
  })
  result1<-result%>%data.frame()%>%t()%>%data.frame()

  netmodx<-cbind(netmod,result1)
  colnames(netmodx)[ncol(netmodx)]<-"thres"

  pp<-ggtexttable(netmodx, theme = ttheme("blank"),rows = NULL) %>%
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
    tab_add_hline(at.row = c(length(netmodx$group)+1), row.side = "bottom", linewidth = 3, linetype = 1) %>%
    tab_add_vline(at.column = 2:(length(colnames(netmodx))), column.side = "left", from.row = 2, linetype = 2) %>%
    tab_add_footnote(text = "*All data processed by R were presented.", size = 10, face = "italic")
  pp

  file2=paste(export_path,"/Network analysis topological roles ( intial version"," ) .txt",sep = "")
  write.table(netmodx,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")

  return(list(net=net,g=g,result=result,netmod=netmodx,plot=pp))
}

"rmnet_matrix_export" <- function(net,export_path="network") {
  export_path<-paste(export_path,"/microbial network analysis/subnet_matrix",sep = "")
  dir.create(export_path, recursive = TRUE)
  for (tt in 1:length(net)) {
    x<-net[[tt]]
    netname<-names(net)[tt]
    cleaned.matrix<-x$cleaned.matrix
    thres<-x$result

    file2=paste(export_path,"/network_matrix (",netname," ",thres,").txt",sep = "")
    write.table(cleaned.matrix,file = file2,  sep = "\t")
    cat("\n","cleaned.matrix (",netname," thres: ",thres," ",nrow(cleaned.matrix)," X ",ncol(cleaned.matrix),")"," has been exported to","/",export_path,"",sep = "","\n")

  }

}

"rmnet_cytogephi_export" <- function(g,taxon,
                                     file.save=TRUE,
                                     export_path ='network/subnet_data_cyto_gephi') {

  taxon$name<-rownames(taxon)
  if (file.save) export_path<-paste(export_path,"/microbial network analysis/subnet_data_cyto_gephi",sep = "")
  if (file.save) dir.create(export_path, recursive = TRUE)


  newdata1<-data.frame()
  newdata3<-data.frame()
  newdata2<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    netname<-names(g)[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    V(g1)$modularity <- membership(cluster_fast_greedy(g1))
    node_table = as.data.frame(vertex.attributes(g1))
    node_table$degree = igraph::degree(g1)
    node_table$closeness=closeness(g1, mode="all")
    node_table$centrality=alpha_centrality(g1,alpha = 0.95)
    node_table$coreness<-coreness(g1)
    node_table$eccentricity<-eccentricity(g1, vids = V(g1), mode = "all")
    node_table$strength<-strength(g1,vids = V(g1),mode = "all",loops = TRUE, weights = NULL)
    node_table$betweenness=betweenness(g1)
    node_table$evcent=evcent(g1)$vector
    module_membership = membership(cluster_fast_greedy(g1))
    node_table$module= membership(cluster_fast_greedy(g1))
    edge_table = as.data.frame(edge.attributes(g1))
    edge_names = unlist(strsplit(attr(E(g1), "vnames"), "|", fixed=T))
    edge_table$node1 = edge_names[seq(from = 1, to = (length(edge_names) - 1), by = 2)]
    edge_table$node2 = edge_names[seq(from = 2, to = length(edge_names), by = 2)]
    edge_table$cor[which(edge_table$correlation < 0)] = "-1"
    edge_table$cor[which(edge_table$correlation > 0)] = "1"
    node_table<-merge(node_table,taxon,by="name")
    node_table$id<-node_table$name
    edge_table$source<-edge_table$node1
    edge_table$target<-edge_table$node2
    edge_table<-subset(edge_table,select = c(6,7,5,2))
    keynode<-zipi(g1)
    keynode[is.na(keynode)]<-0

    keynode$roles<-"Peripheral"
    keynode$roles[which(keynode$pi >= 0.62 & keynode$zi < 2.5)] = "Connector"
    keynode$roles[which(keynode$pi >= 0.62 & keynode$zi >= 2.5)] = "Network hub"
    keynode$roles[which(keynode$pi < 0.62 & keynode$zi >= 2.5)] = "Module hub"

    colnames(keynode)<-c("name","modules","zi","pi",'role')
    node_table<-merge(node_table,keynode,by="name")

    edge<-edge_table
    edge$group<-i
    hub<-node_table
    hub$group<-i
    modu<-V(g1)$modularity
    netcartomod<-max(table(modu))/sum(unique(modu))


if (file.save) {
    file1=paste(export_path,"/", netname, "_cyto_node",".txt",sep = "")
    write.table(node_table,file = file1, row.names = F, quote = F, sep = "\t")

    file1=paste(export_path,"/", netname, "_cyto_edge",".txt",sep = "")
    write.table(edge_table,file = file1, row.names = F, quote = F, sep = "\t")

    flie=paste(export_path,"/", netname, "_gephi_node",".csv", sep = "")
    write.csv(node_table,file = flie,row.names = FALSE)

    flie=paste(export_path,"/", netname, "_gephi_edge",".csv", sep = "")
    write.csv(edge_table,file = flie,row.names = FALSE)
    cat("\n",paste("network",i," (",netname,") ",sep = ""),"used for visualization in cytoscape and gephi has been exported to",export_path,"\n")
}
    {
      newdata1<-rbind(newdata1,netcartomod)
      newdata2<-rbind(newdata2,edge)
      newdata3<-rbind(newdata3,hub)
    }
  }

  newdata3<-newdata3
  newdata2<-newdata2
  return(list(edge_table=newdata2,node_table=newdata3,netname=names(g)))
}

"rmnet_circosfile" <- function(rmnet_cg,abun,
                               export_path ='network/circos') {
  export_path<-paste(export_path,"/microbial network analysis/circos",sep = "")
  dir.create(export_path, recursive = TRUE)
  newdata3<-rmnet_cg$node_table
  newdata2<-rmnet_cg$edge_table
  otu_rare<-abun

  newdata3$group<-as.factor(newdata3$group)

  otusum<-rowSums(otu_rare)%>%data.frame()
  colnames(otusum)<-"value"
  otusum$name <- rownames(otusum)

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
    split_otu <- lapply(apply(sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),2,FUN = function(x){otu_rare[,x]}),function(x){x[-(which(rowSums(x)==0)),]})

  } else {
    split_otu <- lapply(lapply(sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),FUN = function(x){otu_rare[,x]}),function(x){x[-(which(rowSums(x)==0)),]})

  }




  subnet_abun_sum <- lapply(split_otu,function(x){
    subnet_abun_sum<-rowSums(x)%>%data.frame()
    colnames(subnet_abun_sum)<-"value"
    subnet_abun_sum$name <- rownames(subnet_abun_sum)
    return(subnet_abun_sum)
  })

  rmnet_circos_export(subnet_abun_sum,newdata2,newdata3,export_path=export_path)

}

"rmnet_circos_export" <- function(subnet_abun_sum,newdata2,newdata3,export_path) {
  for(k in 1:length(subnet_abun_sum)){
    rmnetname<-names(subnet_abun_sum)[[k]]

    newdata8<-newdata2[newdata2$group==k,]
    newdata8<-subset(newdata8,select = c(1:3))
    colnames(newdata8)<-c("Source","Target","cor")


    newdata4<-newdata3[newdata3$group==k,]
    newdata5<-subnet_abun_sum[[k]]

    newdata6<-merge(newdata4,newdata5,by="name")
    newdata7<-subset(newdata6,select = c(1,12:18,25))

    colnames(newdata7)[1]<-"geneid"
    colnames(newdata7)[9]<-"fc"

    cat("\n","Exporting","",k,paste(" node-circos data.table (",rmnetname,")",sep = ""),"\n")
    cat("\n","Exporting","",k,paste("edge-circos data.table (",rmnetname,")",sep = ""),"\n")

    flie=paste(export_path, "/node-circos",k,".csv", sep = "")
    write.csv(newdata7,file = flie,row.names = FALSE)

    flie=paste(export_path, "/edge-circos",k,".csv", sep = "")
    write.csv(newdata8,file = flie,row.names = FALSE)

  }

}

"rmnet_circosToPerl" <- function(rmnet_cg,taxon,
                                 export_path ='network/circos') {

  options(warn = -1)
  export_path<-paste(export_path,"/microbial network analysis/circos",sep = "")
  dir.create(export_path, recursive = TRUE)
  for (tt in 1:length(unique(rmnet_cg$node_table$group))) {

    select_group<-paste("circos",tt,sep = "")

    df = read.csv(paste(export_path,'/node-',select_group,'.csv',sep = ""),header = T)
    edge_circos<-read.csv(paste(export_path,'/edge-',select_group,'.csv',sep = ""),header = T)
    library(tidyverse)

    df<-circosfile(df)

    ddf<-df[order(df$Kingdom,df$Phylum,
                  df$Class,df$Order,
                  df$Family,df$Genus,
                  df$Species),]



    dir.create(paste(export_path,'/circosabun/',select_group,sep = ""), recursive = TRUE)
    export_path_new=paste(export_path,'/circosabun/',select_group,"/",sep = "")
    circosfile_abun(ddf,taxon,edge_circos,export_path_new)

    cat("\n","Circos plot based on abundance table has been exported to ",export_path_new,"\n")


    dir.create(paste(export_path,'/circoscom/',select_group,sep = ""), recursive = TRUE)
    export_path_new=paste(export_path,'/circoscom/',select_group,"/",sep = "")
    circosfile_compos(ddf,taxon,edge_circos,export_path_new)

    cat("\n","Circos plot based on community structure has been exported to ",export_path_new,"\n")


  }

}

"rmnet_circosColorToPerl" <- function(rmnet_cg,
                 color.circos=colorCustom(50,pal = "gygn"),
                 color.circos.alpha=NULL,
                 export_path ='network/circos') {

  export_path<-paste(export_path,"/microbial network analysis/circos",sep = "")
  dir.create(export_path, recursive = TRUE)
  repnum<-(nrow(rmnet_cg$node_table)/length(color.circos))%>%as.integer()+1
  color.use<-rep(color.circos,repnum)

  rgb<-col2rgb(color.use)%>%t()%>%data.frame()
  colnames(rgb)<-c("rgb1",'rgb2','rgb3')
  if (!is.null(color.circos.alpha)) rgb$alpha<- color.circos.alpha else rgb$alpha<- 1

  rgb$rr<-paste0(rgb$rgb1,",",rgb$rgb2,",",rgb$rgb3,",",rgb$alpha)
  rgb$RGB<-rgb$rr

  rgb$silverx<-paste("silver",1:nrow(rgb),sep = "")
  rgb$silverx<-paste(rgb$silverx,"=",sep = " ")
  silver<-unite(rgb,'col',silverx,RGB,sep = " ")
  silver<-silver$col%>%data.frame()
  write.table(silver, paste(export_path,'/circos_colors.hsv.conf.txt',sep = ""), row.names = FALSE, col.names=FALSE,sep = '\t', quote = FALSE)

  cat("\n","Random colours matching the count of species in the current circos map have been exported to",export_path,"\n")
  message("\n","This file needs to be copied to /Users/cml/software/circos/tutorials/tutorials/10/1/colors.hsv.conf","\n")

}

"rmnet_topodata" <- function(newdata3,abun,taxon,taxanum) {
  otu_rare<-abun
  otusum<-rowSums(otu_rare)%>%data.frame()
  colnames(otusum)<-"value"
  otusum$name <- rownames(otusum)

  phylum<-subset(taxon,select=Phylum)%>%data.frame()
  phylum$name <- rownames(phylum)
  phy <- merge(otusum, phylum, by = "name")
  phy <- subset(phy, select = -name)

  phy %>%
    group_by(Phylum) %>%
    summarise_all(sum) -> phy
  phy <- as.data.frame(phy)
  phy <- na.omit(phy)
  rownames(phy) <- phy$Phylum
  # rowsums and sorted
  phy_table <- phy[order(phy$value,
                         decreasing = TRUE), ]%>%data.frame()
  # top 10
  phy_table10 <- phy_table[1:taxanum, ]%>%data.frame()
  phy_table10['Others', ] <- sum(phy_table$value) - sum(phy_table10$value)


  genus<-subset(taxon,select=Genus)%>%data.frame()
  genus$name <- rownames(genus)
  gen <- merge(otusum, genus, by = "name")
  gen <- subset(gen, select = -name)

  gen %>%
    group_by(Genus) %>%
    summarise_all(sum) -> gen
  gen <- as.data.frame(gen)
  gen <- na.omit(gen)
  rownames(gen) <- gen$Genus

  # rowsums and sorted
  gen_table <- gen[order(gen$value,
                         decreasing = TRUE), ]%>%data.frame()
  # top 10
  gen_table10 <- gen_table[1:taxanum, ]%>%data.frame()
  gen_table10['Others', ] <- sum(gen_table$value) - sum(gen_table10$value)


  data3<-newdata3
  data3$Phylum[!(data3$Phylum %in% phy_table10$Phylum[1:taxanum])]<-"Others"
  data3$Genus[!(data3$Genus %in% gen_table10$Genus[1:taxanum])]<-"Others"


  data3$Phylum<-gsub("p__",'',data3$Phylum)
  data3$Genus<-gsub("g__",'',data3$Genus)
  return(data3)
}

"plot_rmnet_ZPplot" <- function(rmnet_cg,abun,taxon,
                              type=c("mixed","all"),
                              taxa="Phylum",
                              taxanum=10,
                              color_taxa=colorCustom(50,pal = "gygn"),
                              export_path ='network/topological_roles') {

  type<-match.arg(type)
  export_path<-paste(export_path,"/microbial network analysis/topological_roles",sep = "")
  dir.create(export_path, recursive = TRUE)

  newdata3<-rmnet_cg$node_table

  data3 <-rmnet_topodata(newdata3=newdata3,
                         abun=abun,taxon=taxon,
                         taxanum=taxanum)


  colname<-colnames(abun)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)

  for (kk in 1:length(trt_id)) {
    data3$group[which(data3$group==kk)]<-trt_id[kk]
  }

  data3$group<-factor(data3$group,levels = unique(data3$group))
  data4<-filter(data3,data3$zi>2.5 | data3$pi >0.62)

  if (type=="mixed") {
    if ( taxa=="Phylum") {
      color_phylum <- color_taxa[1:length(unique(data3$Phylum))]
      names(color_phylum) <- unique(data3$Phylum)

      library(ggrepel)

      ppp_phy<-ggplot() +
        geom_point(data=data3,mapping=aes(x=pi,y=zi,
                                          size=degree,
                                          color=Phylum,
                                          shape=group,fill=Phylum)) +
        scale_discrete_manual(values=color_phylum,
                              aesthetics = "colour")+
        scale_discrete_manual(values=color_phylum,
                              aesthetics = "fill")+
        annotate(geom = "text",x = -0.1, size = 6,y = 4,
                 label = paste("Module hub","\n","(",nrow(data3[which(data3$role=="Module hub"),]),")"),
                 family = "serif")+
        annotate(geom = "text",x = -0.1, size = 6,y = min(data3$zi),
                 label = paste("Peripheral","\n","(",nrow(data3[which(data3$role=="Peripheral"),]),")"),
                 family = "serif")+
        annotate(geom = "text",x = 0.9, size = 6,y = 4,
                 label = paste("Network hub","\n","(",nrow(data3[which(data3$role=="Network hub"),]),")"),
                 family = "serif")+
        annotate(geom = "text",x = 0.9, size = 6,y = min(data3$zi),
                 label = paste("Connectors","\n","(",nrow(data3[which(data3$role=="Connector"),]),")"),
                 family = "serif")+
        geom_text_repel(data = data4, aes(x = pi, y = zi, label = name),
                        family = "serif",size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)+
        theme_test()+theme(legend.position = 'right')+
        scale_x_continuous(limits = c(-0.2,1),expand = c(0,0))+
        scale_y_continuous(limits = c(min(data3$zi),4))+
        geom_hline(aes(yintercept=2.5),linetype=5,col="red")+
        geom_vline(aes(xintercept=0.62),linetype=5,col="red")+theme(aspect.ratio = 0.8)+
        ylab("Within-module connectivity (Zi)") +
        xlab("Among-module connectivity (Pi)")+
        theme(
          axis.title.x=element_text(colour='black', size=18,face = "bold.italic",
                                    family = "serif", vjust = -1),
          axis.title.y=element_text(colour='black', size=18,face = "bold.italic",
                                    family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=14,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 14,angle = 0,hjust = 0.5,
                                   vjust =0.5,family = "serif"),
          legend.text = element_text(colour='black', size=12,face = "bold",family = "serif"),
          legend.title = element_text(colour='black', size=12,face = "bold",family = "serif"))
      ppp_phy
      ggsave(paste(export_path,'/',"toporole_phylum_level",'.pdf', sep = ''),
             ppp_phy)}

    if ( taxa=="Genus") {
      color_phylum <- color_taxa[1:length(unique(data3$Genus))]
      names(color_phylum) <- unique(data3$Genus)

      library(ggrepel)

      ppp_gen<-ggplot() +
        geom_point(data=data3,mapping=aes(x=pi,y=zi,
                                          size=degree,
                                          color=Genus,
                                          shape=group,fill=Genus)) +
        scale_discrete_manual(values=color_phylum,
                              aesthetics = "colour")+
        scale_discrete_manual(values=color_phylum,
                              aesthetics = "fill")+
        annotate(geom = "text",x = -0.1, size = 6,y = 4,
                 label = paste("Module hub","\n","(",nrow(data3[which(data3$role=="Module hub"),]),")"),
                 family = "serif")+
        annotate(geom = "text",x = -0.1, size = 6,y = min(data3$zi),
                 label = paste("Peripheral","\n","(",nrow(data3[which(data3$role=="Peripheral"),]),")"),
                 family = "serif")+
        annotate(geom = "text",x = 0.9, size = 6,y = 4,
                 label = paste("Network hub","\n","(",nrow(data3[which(data3$role=="Network hub"),]),")"),
                 family = "serif")+
        annotate(geom = "text",x = 0.9, size = 6,y = min(data3$zi),
                 label = paste("Connectors","\n","(",nrow(data3[which(data3$role=="Connector"),]),")"),
                 family = "serif")+geom_text_repel(data = data4, aes(x = pi, y = zi, label = name),
                                                   family = "serif",size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)+
        theme_test()+theme(legend.position = 'right')+
        scale_x_continuous(limits = c(-0.2,1),expand = c(0,0))+
        scale_y_continuous(limits = c(min(data3$zi),4))+
        geom_hline(aes(yintercept=2.5),linetype=5,col="red")+
        geom_vline(aes(xintercept=0.62),linetype=5,col="red")+theme(aspect.ratio = 0.8)+
        ylab("Within-module connectivity (Zi)") +
        xlab("Among-module connectivity (Pi)")+
        theme(
          axis.title.x=element_text(colour='black', size=18,face = "bold.italic",
                                    family = "serif", vjust = -1),
          axis.title.y=element_text(colour='black', size=18,face = "bold.italic",
                                    family = "serif",vjust = 1.5),
          axis.text.y=element_text(colour='black',size=14,family = "serif"),
          axis.text.x=element_text(colour = "black",size = 14,angle = 0,hjust = 0.5,
                                   vjust =0.5,family = "serif"),
          legend.text = element_text(colour='black', size=12,face = "bold",family = "serif"),
          legend.title = element_text(colour='black', size=12,face = "bold",family = "serif"))
      ppp_gen
      ggsave(paste(export_path,'/',"toporole_genus_level",'.pdf', sep = ''),
             ppp_gen)}}
  if (type=="all") {
    if ( taxa=="Phylum") {
      color_phylum <- color_taxa[1:length(unique(data3$Phylum))]
      names(color_phylum) <- unique(data3$Phylum)

      library(ggrepel)

      for (tt in 1:length(trt_id)) {

        ddata3<-data3[which(data3$group==trt_id[tt]),]
        ddata4<-data4[which(data4$group==trt_id[tt]),]

        ppp_phy<-ggplot() +
          geom_point(data=ddata3,mapping=aes(x=pi,y=zi,
                                             size=degree,
                                             color=Phylum,
                                             fill=Phylum)) +
          scale_discrete_manual(values=color_phylum,
                                aesthetics = "colour")+
          scale_discrete_manual(values=color_phylum,
                                aesthetics = "fill")+
          annotate(geom = "text",x = -0.1, size = 6,y = 4,
                   label = paste("Module hub","\n","(",nrow(ddata3[which(ddata3$role=="Module hub"),]),")"),
                   family = "serif")+
          annotate(geom = "text",x = -0.1, size = 6,y = min(ddata3$zi),
                   label = paste("Peripheral","\n","(",nrow(ddata3[which(ddata3$role=="Peripheral"),]),")"),
                   family = "serif")+
          annotate(geom = "text",x = 0.9, size = 6,y = 4,
                   label = paste("Network hub","\n","(",nrow(ddata3[which(ddata3$role=="Network hub"),]),")"),
                   family = "serif")+
          annotate(geom = "text",x = 0.9, size = 6,y = min(ddata3$zi),
                   label = paste("Connectors","\n","(",nrow(ddata3[which(ddata3$role=="Connector"),]),")"),
                   family = "serif")+geom_text_repel(data = ddata4, aes(x = pi, y = zi, label = name),
                                                     family = "serif",size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)+
          theme_test()+theme(legend.position = 'right')+
          scale_x_continuous(limits = c(-0.2,1),expand = c(0,0))+
          scale_y_continuous(limits = c(min(ddata3$zi),4))+
          geom_hline(aes(yintercept=2.5),linetype=5,col="red")+
          geom_vline(aes(xintercept=0.62),linetype=5,col="red")+theme(aspect.ratio = 0.8)+
          ylab("Within-module connectivity (Zi)") +
          xlab("Among-module connectivity (Pi)")+
          theme(
            axis.title.x=element_text(colour='black', size=18,face = "bold.italic",
                                      family = "serif", vjust = -1),
            axis.title.y=element_text(colour='black', size=18,face = "bold.italic",
                                      family = "serif",vjust = 1.5),
            axis.text.y=element_text(colour='black',size=14,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 14,angle = 0,hjust = 0.5,
                                     vjust =0.5,family = "serif"),
            legend.text = element_text(colour='black', size=12,face = "bold",family = "serif"),
            legend.title = element_text(colour='black', size=12,face = "bold",family = "serif"))

        ggsave(paste(export_path,'/',"toporole_phylum_level (",trt_id[tt],').pdf', sep = ''),
               ppp_phy)
      }

    }

    if ( taxa=="Genus") {
      color_phylum <- color_taxa[1:length(unique(data3$Genus))]
      names(color_phylum) <- unique(data3$Genus)

      library(ggrepel)

      for (tt in 1:length(trt_id)) {

        ddata3<-data3[which(data3$group==trt_id[tt]),]
        ddata4<-data4[which(data4$group==trt_id[tt]),]

        ppp_gen<-ggplot() +
          geom_point(data=ddata3,mapping=aes(x=pi,y=zi,
                                             size=degree,
                                             color=Genus,
                                             fill=Genus)) +
          scale_discrete_manual(values=color_phylum,
                                aesthetics = "colour")+
          scale_discrete_manual(values=color_phylum,
                                aesthetics = "fill")+
          annotate(geom = "text",x = -0.1, size = 6,y = 4,
                   label = paste("Module hub","\n","(",nrow(ddata3[which(ddata3$role=="Module hub"),]),")"),
                   family = "serif")+
          annotate(geom = "text",x = -0.1, size = 6,y = min(ddata3$zi),
                   label = paste("Peripheral","\n","(",nrow(ddata3[which(ddata3$role=="Peripheral"),]),")"),
                   family = "serif")+
          annotate(geom = "text",x = 0.9, size = 6,y = 4,
                   label = paste("Network hub","\n","(",nrow(ddata3[which(ddata3$role=="Network hub"),]),")"),
                   family = "serif")+
          annotate(geom = "text",x = 0.9, size = 6,y = min(ddata3$zi),
                   label = paste("Connectors","\n","(",nrow(ddata3[which(ddata3$role=="Connector"),]),")"),
                   family = "serif")+
          geom_text_repel(data = ddata4, aes(x = pi, y = zi, label = name),
                          family = "serif",size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)+
          theme_test()+theme(legend.position = 'right')+
          scale_x_continuous(limits = c(-0.2,1),expand = c(0,0))+
          scale_y_continuous(limits = c(min(ddata3$zi),4))+
          geom_hline(aes(yintercept=2.5),linetype=5,col="red")+
          geom_vline(aes(xintercept=0.62),linetype=5,col="red")+theme(aspect.ratio = 0.8)+
          ylab("Within-module connectivity (Zi)") +
          xlab("Among-module connectivity (Pi)")+
          theme(
            axis.title.x=element_text(colour='black', size=18,face = "bold.italic",
                                      family = "serif", vjust = -1),
            axis.title.y=element_text(colour='black', size=18,face = "bold.italic",
                                      family = "serif",vjust = 1.5),
            axis.text.y=element_text(colour='black',size=14,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 14,angle = 0,hjust = 0.5,
                                     vjust =0.5,family = "serif"),
            legend.text = element_text(colour='black', size=12,face = "bold",family = "serif"),
            legend.title = element_text(colour='black', size=12,face = "bold",family = "serif"))

        ggsave(paste(export_path,'/',"toporole_genus_level (",trt_id[tt],').pdf', sep = ''),
               ppp_gen)
      }
    }
  }


}

"rmnet_keynode_export" <- function(rmnet_cg,
                                   xlabname=c("CS0","CS100","CA100","CD100"),
                                     export_path ='network/topological_roles') {

  export_path<-paste(export_path,"/microbial network analysis/topological_roles",sep = "")
  dir.create(export_path, recursive = TRUE)

  newdata3<-rmnet_cg$node_table
  groupname<-rmnet_cg$netname
if (!is.null(xlabname)) names(groupname)<-xlabname else names(groupname)<-groupname
  for (tk in 1:length(groupname)) {
    newdata3$group[which(newdata3$group==tk)]<-groupname[tk]
  }


  pizi<-filter(newdata3,newdata3$zi>2.5 | newdata3$pi >0.62)

  if (nrow(pizi)==0) {
    stop("\n","No key nodes could be detected !!!")
  }
  tab1<-subset(pizi,select = c(name,group,role,module,degree,Phylum,Class,Order,Family,Genus))
  tab1$group<-factor(tab1$group,levels =rmnet_cg$netname )
  tab1<-tab1[order(tab1$group,tab1$role,tab1$module),]
  tab2<-tab1
  tab1$group<-names(groupname)[match(tab1$group,groupname)]
  tab2$regroup<-names(groupname)[match(tab2$group,groupname)]
  Topological_roles <- ggtexttable(tab1, theme = ttheme("blank"),rows = NULL) %>%
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
    tab_add_hline(at.row = c(length(tab1$name)+1), row.side = "bottom", linewidth = 3, linetype = 1) %>%
    tab_add_vline(at.column = 2:(length(colnames(tab1))), column.side = "left", from.row = 2, linetype = 2) %>%
    tab_add_footnote(text = "*All data processed by R were presented.", size = 10, face = "italic")

  #ggsave('~/Desktop/11/all/topological_roles/Topological_roles.pdf')
  print(Topological_roles)
  flie=paste(export_path,"/keynode",".csv", sep = "")
  write.csv(tab1,file = flie,row.names = FALSE)
  cat("\n","Key nodes topological roles result has been exported to ",export_path,"\n")

  return(tab2)
}

"powerlaw_fit" <- function(g1,decimal=3,exchantime=999) {

  V(g1)$degree <- igraph::degree(g1)

  degree_dist <- table(V(g1)$degree)
  degree_num <- as.numeric(names(degree_dist))
  degree_count <- as.numeric(degree_dist)

  dat <- data.frame(degree = degree_num, count = degree_count)
  intial<-max(dat$count)
  lambda<-(dat$count[3]-intial)/dat$count[3]
  mod <- stats::nls(dat$count ~ a*dat$degree^b, data = dat,
                    algorithm="default",
                    start = list(a = intial, b = lambda))


  summary(mod)

  fit <- fitted(mod)
  SSre <- sum((dat$count-fit)^2)
  SStot <- sum((dat$count-mean(dat$count))^2)
  R2 <- round(1 - SSre/SStot, decimal)
  R2

  a <- round(coef(mod)[1], decimal)
  b <- round(coef(mod)[2], decimal)
  a; b

  p_num <- 1
  dat_rand <- dat
  for (i in 1:exchantime) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count-fit)^2)
    SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
  }
  p_value <- p_num / (exchantime+1)
  p_value

  label <- data.frame(
    formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
    R2 = sprintf('italic(R^2) == %.3f', R2),
    p_value = sprintf('italic(P) < %.3f', p_value)
  )



  return(list(a=a,lambda=b,r2=R2,p=p_value,formula=label))
}

"powerlaw_logfit" <- function(g1,decimal=3,exchantime=999) {

  V(g1)$degree <- igraph::degree(g1)

  degree_dist <- table(V(g1)$degree)
  degree_num <- as.numeric(names(degree_dist))
  degree_count <- as.numeric(degree_dist)

  dat <- data.frame(degree = log(degree_num), count = log(degree_count))

  intial<-max(dat$count)
  lambda<-(dat$count[4]-intial)/dat$degree[4]

  mod <- nls(dat$count ~ b*dat$degree+a, data = dat,
             algorithm="default",
             start = list(a = intial,
                          b = lambda))


  summary(mod)

  fit <- fitted(mod)
  SSre <- sum((dat$count-fit)^2)
  SStot <- sum((dat$count-mean(dat$count))^2)
  R2 <- round(1 - SSre/SStot, decimal)
  R2
  a <- round(coef(mod)[1], decimal)
  b <- round(coef(mod)[2], decimal)
  a; b

   p_num <- 1
  dat_rand <- dat
  for (i in 1:exchantime) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count-fit)^2)
    SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
  }
  p_value <- p_num / (exchantime+1)
  p_value

  label <- data.frame(
    formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
    R2 = sprintf('italic(R^2) == %.3f', R2),
    p_value = sprintf('italic(P) < %.3f', p_value)
  )


  return(list(a=a,lambda=b,r2=R2,p=p_value,formula=label))
}

"randnet" <- function(g1,times=1000,decimal=3) {

  node_num <- vcount(g1)
  edge_num <- ecount(g1)
  V(g1)$degree <- igraph::degree(g1)
  average_degree <- mean(igraph::degree(g1))
  clustering_coefficient <- transitivity(g1,type="average", isolates="zero")    #聚类系数
  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)
  modularity<-modularity(cluster_fast_greedy(g1))
  dist <- mean_distance(g1, directed=F, unconnected=T)

  degree_dist <- table(V(g1)$degree)
  degree_num <- as.numeric(names(degree_dist))
  degree_count <- as.numeric(degree_dist)
  names(degree_count) <- degree_num
  degs <- rep(degree_num, degree_count)

  clustering_coefficient_rand <- clustering_coefficient
  modularity_rand<-modularity
  dist_rand <- dist
  average_degree_rand <- average_degree
  result <- NULL
  set.seed(123)
  for (i in 1:times) {
    g_rand <- degree.sequence.game(degs, method = 'vl')
    clustering_coefficient_rand <- c(clustering_coefficient_rand,
                                     transitivity(g_rand))
    modularity_rand<-c(modularity_rand,modularity(walktrap.community(g_rand)))

    dist_rand<-c(dist_rand,mean_distance(g_rand, directed=F, unconnected=T))

  }

  result <- rbind(result, c(mean(clustering_coefficient_rand),
                            sd(clustering_coefficient_rand),
                            sd(clustering_coefficient_rand)/sqrt(times)),
                  c(mean(modularity_rand),
                    sd(modularity_rand),
                    sd(modularity_rand)/sqrt(times)),
                  c(mean(dist_rand),
                    sd(dist_rand),
                    sd(dist_rand)/sqrt(times))

  )

  result <- data.frame(result)
  colnames(result) <- c("mean","sd","se")
  rownames(result) <- c("rand average cluster coefficient (avgCC)",
                        "rand modularity (M)",
                        "rand average path distance (GD)")
  decimal
  result$res<-paste(sprintf("%0.3f",round(result$mean, 3)),
                    sprintf("%0.3f",round(result$se, 3)),
                    sep = '±')

  res<-t(result$res)%>%data.frame()
  colnames(res)<-c("avgCC",
                   "modularity",
                   "GD")


  cat("\n","Calculating",times," random networks and their topological roles","\n")
  return(res)
}

"rmnet_toporoles_export" <- function(rmnet_igraph,randcalc_times=100,
                                     add.params=NULL,
                                     xlabname=c("CS0","CS100","CA100","CD100"),
                                     export_path ='network/topological_roles') {
  export_path<-paste(export_path,"/microbial network analysis/topological_roles",sep = "")
  dir.create(export_path, recursive = TRUE)
  g<-rmnet_igraph$g
  result<-rmnet_igraph$result

  newdata<-data.frame()

  for(i in 1:length(g)){

    g1 <- g[[i]]

    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    gd = cluster_fast_greedy(g1) # greedy module separation
    gd
    modules = length(gd) # number of modules
    modules
    M = modularity(gd) # modularity
    M
    avgK = mean(centr_degree(g1)$res) # average degree
    avgK
    avgCC = transitivity(g1, type="average", isolates="zero") # average clustering coefficient. Zero for bipartite network (no triangular subnetwork).
    avgCC
    GD = mean_distance(g1, directed=F, unconnected=T) # geodesic distance
    GD


    edge_table = as.data.frame(edge.attributes(g1))
    edge_table$cor[which(edge_table$correlation < 0)] = "-1"
    edge_table$cor[which(edge_table$correlation > 0)] = "1"
    neg<-length(edge_table$cor[which(edge_table$cor == -1)])/length(edge_table$cor)*100
    pos<-length(edge_table$cor[which(edge_table$cor == 1)])/length(edge_table$cor)*100

    #powerfit<-powerlaw_fit(g1)
    #powerlaw_logfit<-powerlaw_logfit(g1)

    #powlaw_r2<-powerfit$r2
    #powlaw_lambda<-powerfit$lambda*(-1)
   # powlaw_p<-powerfit$p
   # powerlaw<-paste("r2: ",powlaw_r2," (lambda: ",powlaw_lambda,", p: ",powlaw_p,")")

    #powlaw_log_r2<-powerlaw_logfit$r2
    #powlaw_log_lambda<-powerlaw_logfit$lambda*-1
    #powlaw_log_p<-powerlaw_logfit$p
   # powerlaw_log<-paste(powlaw_log_r2," (",powlaw_log_lambda,", ",powlaw_log_p,")")


    re<-result[[i]]

    {
      topo<-rbind(
        re,
        #powerlaw,
        #powerlaw_log,
        length(igraph::degree(g1)),
        nrow(data.frame(as_edgelist(g1))),
        sprintf("%0.2f",round(pos,2)),
        sprintf("%0.2f",round(neg,2)),
        sprintf("%0.3f",round(avgK,3)),
        modules,
        sprintf("%0.3f",round(avgCC,3)),
        sprintf("%0.3f",round(M,3)),
        sprintf("%0.3f",round(GD,3))
      )%>%data.frame()
      topo<-t(topo)%>%data.frame()
      newdata<-rbind(newdata,topo)
    }
  }
  newdata<-t(newdata)%>%data.frame()
  newdata$name<-c('Similarity threshold (St)',
                  #'powerlaw',
                 # 'powerlaw_log',
                  'Node number (N)',
                  'Edge number (N)',
                  'Positive edges (%)',
                  'Negative edges (%)',
                  'Average degree (avgK)',
                  'Module number (N)',
                  'Average cluster coefficient (avgCC)',
                  'Modularity (M)',
                 'Average path distance (GD)')
  newdata<-subset(newdata,select = c(name,1:length(g)))
  colnames(newdata)<-c('parameters',xlabname)

  newdata2<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    dft<-randnet(g1,times = randcalc_times)%>% data.frame()

    {
      newdata2<-rbind(newdata2,dft)
    }
  }


  newdata2<-t(newdata2)%>%data.frame()
  newdata2$parameters<-rownames(newdata2)
  newdata2<-subset(newdata2,select = c(parameters,1:length(g)))
  colnames(newdata2)<-c('parameters',xlabname)
  newdata2$parameters<-c('Average cluster coefficient (avgCC)',
                         'Modularity (M)',
                         'Average path distance (GD)')
  netn<-data.frame(rep("",1+length(g)))%>%t()%>%data.frame()
  netn$X1[1]<-"Empirical Network"
  colnames(netn)<-colnames(newdata)
  newdata<-rbind(netn,newdata)

  if (!is.null(add.params)) for (kk in 1:length(add.params)) {
    if (class(add.params[[kk]])=="data.frame") {
      add.params[[kk]]<-add.params[[kk]]$value[match(add.params[[kk]]$group,names(g))]
    }
  }


  if (!is.null(add.params)) {
    para.sel<-data.frame()
    for (jk in 1:length(add.params)) {
    par.sel<-sprintf("%0.3f",round(add.params[[jk]],3))
    par.sel<-c(names(add.params)[jk],par.sel)
    {
      para.sel<-rbind(para.sel,par.sel)
    }
    }
    colnames(para.sel)<-colnames(newdata)
    newdata<-rbind(newdata[1:8,],para.sel,newdata[9:nrow(newdata),])
  }

  netr<-netn
  netr$parameters[1]<-"Random Network"
  newdata2<-rbind(netr,newdata2)
  newdata<-rbind(newdata,newdata2)
rownames(newdata)<-1:nrow(newdata)



  emnet_randnet_com <- ggtexttable(newdata, theme = ttheme("blank"),rows = NULL) %>%
    tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
    tab_add_hline(at.row = c(length(newdata$parameters)+1), row.side = "bottom", linewidth = 3, linetype = 1) %>%
    tab_add_vline(at.column = 2:(length(g)+1), column.side = "left", from.row = 2, linetype = 2) %>%
    tab_add_footnote(text = "*All data processed by R were presented.", size = 10, face = "italic")

  flie=paste(export_path,"/Topological_roles",".csv", sep = "")
  write.csv(newdata,file = flie,row.names = FALSE)
  cat("\n","All co-occurrence netoworks topological roles result have been exported to ",export_path,"下")
print(emnet_randnet_com)
  return(newdata)
}

"nodechar" <- function(g1,boot=1000) {
  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)

  all_degree <- igraph::degree(g1, mode = 'all')
  all_transitivity <- transitivity(g1, 'local', vids = V(g1))
  all_transitivity[is.na(all_transitivity)] <- 0
  names(all_transitivity)<- V(g1)$name
  all_eccentricity<-eccentricity(g1, vids = V(g1), mode = "all")
  all_betweenness<-betweenness(g1)

  set.seed(123)
  boot_degree <- replicate(boot, sample(all_degree, 1, replace = TRUE))
  boot_transitivity <- replicate(boot, sample(all_transitivity, 1, replace = TRUE))
  boot_eccentricity <- replicate(boot, sample(all_eccentricity, 1, replace = TRUE))
  boot_betweenness <- replicate(boot, sample(all_betweenness, 1, replace = TRUE))

  node_characterstics <- data.frame(cbind(boot_degree,
                                          boot_transitivity,
                                          boot_eccentricity,
                                          boot_betweenness))
  return(node_characterstics)
}

"scale_free_esti" <- function(network_node) {
  length(network_node)

  par(mfrow = c(3, length(network_node)))
  for (network in 1:length(network_node)) {
    hist(network_node[[network]][ ,'boot_degree'], xlab = 'degree',
         ylab = 'count', main = paste0('Degree distribution'," (", network,")"))
  }

  for (network in 1:length(network_node)) {
    degree_dist <- table(network_node[[network]][ ,'boot_degree'])
    degree_num <- as.numeric(names(degree_dist))
    degree_count <- as.numeric(degree_dist)
    plot(degree_num, degree_count, xlab = 'Degree', ylab = 'Count',
         main = paste0('Degree distribution'," (", network,")"))
  }
  for (network in 1:length(network_node)) {
    degree_dist <- table(network_node[[network]][ ,'boot_degree'])
    degree_num <- as.numeric(names(degree_dist))
    degree_count <- as.numeric(degree_dist)
    plot(log(degree_num), log(degree_count), xlab = 'log-degree',
         ylab = 'log-count', main = paste0('Degree distribution'," (", network,")"))

  }

}

"nettopo_com_kstest" <- function(network_node) {

  result <- NULL
  num <- length(network_node)

  for (i in 1:(num-1)) {
    for (j in (i+1):num) {

      node_degree_i <- network_node[[1]][ ,'boot_degree']
      node_degree_j <- network_node[[2]][ ,'boot_degree']
      ks_degree <- ks.test(node_degree_i, node_degree_j, alternative = 'two.sided')


      node_betweenness_i <- network_node[[i]][ ,'boot_betweenness']
      node_betweenness_j <- network_node[[j]][ ,'boot_betweenness']
      ks_betweenness <- ks.test(node_betweenness_i, node_betweenness_j, alternative = 'two.sided')


      node_eccentricity_i <- network_node[[i]][ ,'boot_eccentricity']
      node_eccentricity_j <- network_node[[j]][ ,'boot_eccentricity']
      ks_eccentricity <- ks.test(node_eccentricity_i, node_eccentricity_j, alternative = 'two.sided')


      node_transitivity_i <- network_node[[i]][ ,'boot_transitivity']
      node_transitivity_j <- network_node[[j]][ ,'boot_transitivity']
      ks_transitivity <- ks.test(node_transitivity_i, node_transitivity_j, alternative = 'two.sided')

      result <- rbind(result, c('group1' = i, 'group2' = j,
                                'ks_degree' = ks_degree$p.value, 'ks_betweenness' = ks_betweenness$p.value,
                                'ks_eccentricity' = ks_eccentricity$p.value, 'ks_transitivity' = ks_transitivity$p.value))
    }
  }

  result <- data.frame(result, stringsAsFactors = FALSE)
  result$ks_degree <- as.numeric(result$ks_degree)
  result$ks_betweenness <- as.numeric(result$ks_betweenness)
  result$ks_eccentricity <- as.numeric(result$ks_eccentricity)
  result$ks_transitivity <- as.numeric(result$ks_transitivity)
  result

  cat("\n","ks.test referred to all topological roles have been calculated.","\n")
  return(result)
}

"rmnet_topo_compare" <- function(g,boot_times=999,
                                 export_path ='network/topo_esti_compare') {
  export_path<-paste(export_path,"/microbial network analysis/topo_esti_compare",sep = "")
  dir.create(export_path, recursive = TRUE)
  ##########

  network_node <- list()
  for(i in 1:length(g)){
    network_node[i]<-list(nodechar(g[[i]],boot = boot_times))
  }

  scale_free_esti(network_node)

  pdf(paste0(export_path,"/boot_times(",boot_times,") scale_free_esti.pdf"), encoding="MacRoman", width=11, height=9)
  scale_free_esti(network_node)
  dev.off()


  net_com<-nettopo_com_kstest(network_node)

  flie=paste(export_path, "/toporoles_compare",".csv", sep = "")
  write.csv(net_com,file = flie,row.names = FALSE)

  cat("\n","subsample nodes about",boot_times,"times and export the chart of node-and-number","\n")
}
