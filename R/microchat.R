"setMicrochat" <- function(otu_table,taxon_table,tree=NULL,params = NULL) {

  suppressMessages(library(ape))
  suppressMessages(library(tidyverse))
  suppressMessages(library(tidytree))
  samplenum<-length(unique(substr(colnames(otu_table),start=1,stop=2)))

  taxon_table<-taxonadj(taxon_table)

  "is.letter" <- function(char) {

    if (!is.na(match(tolower(char), letters))) {
      print(TRUE)
    } else {
      print(FALSE)
    }

  }

  if (is.letter(substr(colnames(otu_table),start = 3,stop = 3)[1])) {
    stop("\n","Sample names must be in the form of the mixture of the group names (two letters) and biological replicates number, such as AA1...AA12...bk1...bk12.")
  }

  if (nrow(otu_table)!=nrow(taxon_table)) {
    message("\n","otu_table didn't shared the same OTUs number with taxon_table")
  }


 if (!is.null(params)) {
  if (class(params)!="list") {
   if (nrow(params)!=ncol(otu_table)) {
     stop("\n","Params_table that provided didn't shared the same sample names with otu_table")
   }
  }

   microchatobj <- list(otu_table,taxon_table,tree,params)

   class(microchatobj) <- "microchat"

   names(microchatobj)[[1]]<-"otu_table"
   names(microchatobj)[[2]]<-"taxon_table"
   names(microchatobj)[[3]]<-"tree"
   names(microchatobj)[[4]]<-"param_table"

   cat("\n","otu_table: ",ncol(microchatobj[[1]]), " samples with ",samplenum," groups", sep = "","\n" )
   cat("\n","taxon_table: ",nrow(microchatobj[[2]]), " OTUs with ", ncol(microchatobj[[2]]), " taxa ", sep = "","\n" )
   if (class(params)!="list") cat("\n","param_table: ",nrow(microchatobj[[4]]), " samples with ", ncol(microchatobj[[4]]), " parameters ", sep = "","\n" )
   if (class(params)=="list") {
     cat("\n", length(microchatobj[[4]])," param_table(s): \t", sep = "")
     for (tt in 1:length(microchatobj[[4]])) {
       cat("\n  ",names(microchatobj[[4]])[tt],": ",
           nrow(microchatobj[[4]][[tt]]), " samples with ", ncol(microchatobj[[4]][[tt]]), " parameters ",
           "\t",sep = "")
     }
   }

   print(microchatobj[[3]])
 } else {
   microchatobj <- list(otu_table,taxon_table,tree,params)

   class(microchatobj) <- "microchat"

   names(microchatobj)[[1]]<-"otu_table"
   names(microchatobj)[[2]]<-"taxon_table"
   names(microchatobj)[[3]]<-"tree"
   names(microchatobj)[[4]]<-"param_table"

   cat("\n","otu_table: ",ncol(microchatobj[[1]]), " Samples with ",samplenum," groups", sep = "","\n" )
   cat("\n","taxon_table: ",nrow(microchatobj[[2]]), " OTUs with ", ncol(microchatobj[[2]]), " taxa ", sep = "","\n" )
   cat("\n","param_table: ", " Not provided. ", sep = "","\n" )

   print(microchatobj[[3]])
 }


  return(microchatobj)
}

"setParamchat" <- function(otu_table=NULL,taxon_table=NULL,tree=NULL,params = NULL) {
  if (!is.null(params) && is.null(otu_table) && is.null(taxon_table) && is.null(tree)) {
    data(abun)
    data(taxon)
    abunxx<-abun

    if (class(params) =="list") paramsx<-params[[1]]
    if (class(params)!="list") paramsx<-params

    for (tkt in 1:100) {
      while(nrow(paramsx)>ncol(abunxx)) {
        abunxx<-cbind(abunxx,abun)
      }
    }
    abunx<-abunxx[1:10,1:nrow(paramsx)]
    colnames(abunx)<-rownames(paramsx)
    taxon<-taxon[1:nrow(abunx),]
    abun<-abunx
    otu_table<-abun
    taxon_table<-taxon
  }

  "is.letter" <- function(char) {

    if (!is.na(match(tolower(char), letters))) {
      print(TRUE)
    } else {
      print(FALSE)
    }

  }


  samplenum<-length(unique(substr(colnames(otu_table),start=1,stop=2)))
  if (is.letter(substr(colnames(otu_table),start = 3,stop = 3)[1])) {
    stop("\n","Sample names must be in the form of the mixture of the group names (two letters) and biological replicates number, such as AA1...AA12...bk1...bk12.")
  }
  taxon_table<-taxonadj(taxon_table)
  if (nrow(otu_table)!=nrow(taxon_table)) {
    message("\n","otu_table didn't shared the same OTUs number with taxon_table")
  }

  if (!is.null(params)) {
    if (class(params)!="list") {
      if (!is.null(otu_table) && !is.null(taxon_table) && nrow(params)!=ncol(otu_table)) {
        stop("\n","Params_table that provided didn't shared the same sample names with otu_table")
      }
    }

    microchatobj <- list(otu_table,taxon_table,tree,params)

    class(microchatobj) <- "microchat"

    names(microchatobj)[[1]]<-"otu_table"
    names(microchatobj)[[2]]<-"taxon_table"
    names(microchatobj)[[3]]<-"tree"
    names(microchatobj)[[4]]<-"param_table"

    cat("\n","otu_table: ",ncol(microchatobj[[1]]), " samples with ",samplenum," groups", sep = "","\n" )
    cat("\n","taxon_table: ",nrow(microchatobj[[2]]), " OTUs with ", ncol(microchatobj[[2]]), " taxa ", sep = "","\n" )
    if (class(params)!="list") cat("\n","param_table: ",nrow(microchatobj[[4]]), " samples with ", ncol(microchatobj[[4]]), " parameters ", sep = "","\n" )
    if (class(params)=="list") {
      cat("\n", length(microchatobj[[4]])," param_table(s): \t", sep = "")
      for (tt in 1:length(microchatobj[[4]])) {
        cat("\n  ",names(microchatobj[[4]])[tt],": ",
            nrow(microchatobj[[4]][[tt]]), " samples with ", ncol(microchatobj[[4]][[tt]]), " parameters ",
            "\t",sep = "")
      }
    }

    print(microchatobj[[3]])
  } else {
    microchatobj <- list(otu_table,taxon_table,tree,params)

    class(microchatobj) <- "microchat"

    names(microchatobj)[[1]]<-"otu_table"
    names(microchatobj)[[2]]<-"taxon_table"
    names(microchatobj)[[3]]<-"tree"
    names(microchatobj)[[4]]<-"param_table"

    cat("\n","otu_table: ",ncol(microchatobj[[1]]), " Samples with ",samplenum," groups", sep = "","\n" )
    cat("\n","taxon_table: ",nrow(microchatobj[[2]]), " OTUs with ", ncol(microchatobj[[2]]), " taxa ", sep = "","\n" )
    cat("\n","param_table: ", " Not provided. ", sep = "","\n" )

    print(microchatobj[[3]])
  }


  return(microchatobj)
}

"taxonadj" <- function(taxon) {

  taxon_table<-taxon
  taxon_table[is.na(taxon_table)]<-""
  classnum<-length(colnames(taxon_table))

  if (classnum==7) {
    colnames(taxon_table)<-tolower(colnames(taxon_table))
    match.colnames<-c("kingdom","phylum","class","order","family","genus","species")

    if(length(setdiff(colnames(taxon_table),match.colnames))>0) colnames(taxon_table)<-match.colnames
    colnames(taxon_table)<-stringr::str_to_title(colnames(taxon_table))
  } else {
    stop("\n","Pleased provide senven classifications (From Kingdom To Species) !!!") }

  charclass<-taxon_table$Phylum%>%unique()
  charclass<-charclass[2]

  if (substr(charclass,start = 2,stop = 3)!="__") {
    for (tt in 1:length(colnames(taxon_table))) {
      for (jj in 1:nrow(taxon_table)) {
        if (taxon_table[jj,tt]=="") {
          addtext<-paste(tolower(substr(colnames(taxon_table)[tt],start = 1,stop = 1)),"__",sep = "")
          taxon_table[jj,tt]<-paste(addtext,sep = "")
        }
      }
      addtext<-paste(tolower(substr(colnames(taxon_table)[tt],start = 1,stop = 1)),"__",sep = "")
      taxon_table[,tt]<-paste(addtext,taxon_table[,tt],sep = "")
    }
  } else {
    for (tt in 1:length(colnames(taxon_table))) {
      for (jj in 1:nrow(taxon_table)) {
        if (taxon_table[jj,tt]=="") {
          addtext<-paste(tolower(substr(colnames(taxon_table)[tt],start = 1,stop = 1)),"__",sep = "")
          taxon_table[jj,tt]<-paste(addtext,sep = "")
        }
      }

    }
    taxon_table<-taxon_table
  }

  return(taxon_table)
}

"taxon_filter" <- function(tree.use,taxon,abun, select_taxa=FALSE) {
  abun$sum<-rowSums(abun)
  abun<-abun[which(abun$sum != 0),]

  taxon$name<-rownames(taxon)
  abun$label<-rownames(abun)
  if (!select_taxa) {
    taxon<-dplyr::right_join(taxon,abun,by=c('name'='label'))
    taxon<-subset(taxon,select=c(1:8))
  }
  if (select_taxa) {
    taxon<-dplyr::left_join(taxon,abun,by=c('name'='label'))
    abun <-subset(taxon,select=-(1:7))
    rownames(abun)<-abun$name
    abun<-subset(abun,select = -c(name, sum))
    abun[is.na(abun)]<-0
    abun$sum<-rowSums(abun)
    abun<-abun[which(abun$sum != 0),]

    taxon<-subset(taxon,select=c(1:8))
    abun$label<-rownames(abun)
    taxon<-dplyr::right_join(taxon,abun,by=c('name'='label'))
    taxon<-subset(taxon,select=c(1:8))
  }


  taxon$Phylum[taxon$Phylum == "" ]<-"p__unclassified"
  phy<-unique(taxon$Phylum)

  group_info<-split(taxon$name,taxon$Phylum)

 if (!is.null(tree.use)) {
   tree11<-microchat::groupOTU.phylo(tree.use,group_info)
  data_tree <-tree11%>%as_tibble()

  to_drop <- data_tree[which(data_tree$group == 0),]$label
  tree11 <- treeio::drop.tip(tree11, to_drop)


  group_info<-split(taxon$name,taxon$Phylum)
  tree11<-microchat::groupOTU.phylo(tree11,group_info)
  data_tree <-tree11%>%as_tibble()

  data_tree<- data_tree[which(data_tree$group != 0),]
  } else {
    tree11<-NULL
    data_tree<-NULL
  }

  group_info1<-reshape2::melt(group_info)
  colnames(group_info1)<-c("name","type")
  group_info2<-merge(taxon,group_info1,by="name")

  abun<-subset(abun,select = c(-label,-sum))
  rownames(taxon)<-taxon$name
  taxon<-subset(taxon,select = c(-name))
  return(list(data_tree=data_tree,tree=tree11,group_info2=group_info2,taxon=taxon,abun=abun))
}

"tidyMicrochat" <- function(mchatobj,
                            group.keep=NULL,
                            group.rm=NULL,
                            sample.keep=NULL,
                            sample.rm=NULL,
                            taxon.keep=NULL,
                            taxon.rm=NULL,
                            filter_pollution=TRUE,
                            filter_Archaea=TRUE) {

  if (class(mchatobj)!="microchat") {
    stop("\n","Please use function 'setMicrochat' to convert the data into a 'microchat' object
")
  }
  otu_table <- mchatobj$otu_table
  taxon_table <- mchatobj$taxon_table
  tree <- mchatobj$tree
  param_table<-mchatobj$param_table

  group_infoxx<-split(rownames(taxon_table),taxon_table$Phylum)

  #kkk<-groupOTU.phylo(tree, group_infoxx)
  ###1.remove pollution: Chloroplast & Mitochondria
  if (filter_pollution) {
    taxon_table <- taxon_table[which( taxon_table$Family != "f__Mitochondria" & taxon_table$Order != "o__Chloroplast" & taxon_table$Class != "c__Chloroplast"),]
  } else {
    taxon_table <- taxon_table
  }

  ###2.remove Archaea
  if (filter_Archaea) {
    taxon_table<-taxon_table[which(taxon_table$Kingdom!="k__Archaea"),]
  } else {
    taxon_table <- taxon_table
  }

###3.remove null phylum
    taxon_table<-taxon_table[which(taxon_table$Phylum!="p__"),]

  if (is.null(group.keep)) {
    taxon_table<-taxon_table
  } else {
    otu_table <- otu_table[,which(substr(colnames(otu_table),start = 1,stop = 2) %in% group.keep)]
    if (class(param_table)=="list") {
      for (tt in 1:length(param_table)) {
        param_table[[tt]] <- param_table[[tt]][which(substr(rownames(param_table[[tt]]),start = 1,stop = 2) %in% group.keep),]
      }
    } else {
      param_table <- param_table[which(substr(rownames(param_table),start = 1,stop = 2) %in% group.keep),]
    }
    cat("\n","Samples belonging to '",group.keep,"' have (has) been kept in microchat object.",sep=" ")
  }

  if (is.null(group.rm)) {
    taxon_table<-taxon_table
  } else {
    otu_table <- otu_table[,-which(substr(colnames(otu_table),start = 1,stop = 2) %in% group.rm)]
    if (class(param_table)=="list") {
      for (tt in 1:length(param_table)) {
        param_table[[tt]] <- param_table[[tt]][-which(substr(rownames(param_table[[tt]]),start = 1,stop = 2) %in% group.rm),]
      }
    } else {
      param_table <- param_table[-which(substr(rownames(param_table),start = 1,stop = 2) %in% group.rm),]
    }
    cat("\n","Samples belonging to '",group.rm,"' have (has) been removed in microchat object.",sep=" ")
  }

  if (is.null(sample.keep)) {
    taxon_table<-taxon_table
  } else {
    otu_table <- otu_table[,which(colnames(otu_table) %in% sample.keep)]
    if (class(param_table)=="list") {
      for (tt in 1:length(param_table)) {
        param_table[[tt]] <- param_table[[tt]][which(rownames(param_table[[tt]]) %in% sample.keep),]
      }
    } else {
      param_table <- param_table[which(rownames(param_table) %in% sample.keep),]
    }
    cat("\n","Samples belonging to '",sample.keep,"' have (has) been kept in microchat object.",sep=" ")
  }

  if (is.null(sample.rm)) {
    taxon_table<-taxon_table
  } else {
    otu_table <- otu_table[,-which(colnames(otu_table) %in% sample.rm)]
    if (class(param_table)=="list") {
      for (tt in 1:length(param_table)) {
        param_table[[tt]] <- param_table[[tt]][-which(rownames(param_table[[tt]]) %in% sample.rm),]
      }
    } else {
      param_table <- param_table[-which(rownames(param_table) %in% sample.rm),]
    }
    cat("\n","Samples belonging to '",sample.rm,"' have (has) been removed in microchat object.",sep=" ")
  }

  allchar<-character()
  for (tk in 1:length(colnames(taxon_table))) {
    allcharr<-taxon_table[,tk]%>%unique()
    {
      allchar<-c(allchar,allcharr)
    }
  }

  if (is.null(taxon.keep)) {
    taxon_table<-taxon_table
  } else {
    taxon.keep1<-taxon.keep[1]
    if (substr(taxon.keep1,start = 2,stop = 3)!="__") {
      taxon.keep<-stringr::str_to_title(taxon.keep)
    } else {
      splitdata<-str_split_fixed(taxon.keep,pattern="__",2)
      splitdata[,2]<-stringr::str_to_title(splitdata[,2])
      taxon.keep<-paste(splitdata[,1],"__",splitdata[,2],sep = "")
    }

    fulltaxon<-lapply(
      lapply(
        sapply(as.list(taxon.keep), grep,allchar),
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

    cat("\n","Taxa belonging to '",fulltaxon,"' have (has) been kept in microchat object.",sep=" ")
  }

  if (is.null(taxon.rm)) {
    taxon_table<-taxon_table
  } else {
    taxon.rm1<-taxon.rm[1]
    if (substr(taxon.rm1,start = 2,stop = 3)!="__") {
      taxon.rm<-stringr::str_to_title(taxon.rm)
    } else {
      splitdata<-str_split_fixed(taxon.rm,pattern="__",2)
      splitdata[,2]<-stringr::str_to_title(splitdata[,2])
      taxon.rm<-paste(splitdata[,1],"__",splitdata[,2],sep = "")
    }

    fulltaxon<-lapply(
      lapply(
        sapply(as.list(taxon.rm), grep,allchar),
        function (x) {
          x<-x[1]
        }),
      function(y) {
        allchar[y]
      })%>%as.character()

    removepos<-NULL
    for (sk in 1:length(colnames(taxon_table))) {
      removepos1<-which(taxon_table[,sk] %in% fulltaxon)
      {
        removepos<-c(removepos,removepos1)
      }
    }

    taxon_table<-taxon_table[-removepos,]

    cat("\n","Taxa belonging to '",fulltaxon,"' have (has) been removed in microchat object.",sep=" ")
  }

  name_select<-rownames(taxon_table)
  otu_table <- otu_table[which(rownames(otu_table) %in% name_select),]

  tree_filter<-taxon_filter(tree.use=tree,taxon=taxon_table,abun=otu_table)
  otu_table<-tree_filter$abun
  taxon_table<-tree_filter$taxon
  tree<-tree_filter$tree

  tidymicrochatobj <- list(otu_table,taxon_table,tree,param_table)

  class(tidymicrochatobj) <- "microchat"

  names(tidymicrochatobj)[[1]]<-"otu_table"
  names(tidymicrochatobj)[[2]]<-"taxon_table"
  names(tidymicrochatobj)[[3]]<-"tree"
  names(tidymicrochatobj)[[4]]<-"param_table"

  samplenum<-length(unique(substr(colnames(otu_table),start=1,stop=2)))

  cat("\n","otu_table: ",ncol(tidymicrochatobj[[1]]), " Samples with ",samplenum," groups", sep = "","\n" )
  cat("\n","taxon_table: ",nrow(tidymicrochatobj[[2]]), " OTUs with ", ncol(tidymicrochatobj[[2]]), " taxa ", sep = "","\n" )

  print(tidymicrochatobj[[3]])
  return(tidymicrochatobj)
}


"sub2sample" <- function(otudata,taxon_table,sample.size=NULL) {
  suppressMessages(library(phyloseq))
  options(warn = -1)
  set.seed(8084)
  otu1 = phyloseq::otu_table(otudata, taxa_are_rows = T)
  phyloseq = phyloseq::phyloseq(otu1)

  if (!is.null(sample.size)) {
    minsize<-sample.size
  } else {
    minsize<-phyloseq::sample_sums(phyloseq)%>%min()
  }


  rare.data = phyloseq::rarefy_even_depth(phyloseq,replace = TRUE, sample.size=minsize)

  rare.otu = rare.data@.Data %>% as.data.frame()
  otu_table<-rare.otu

  rare.otu$name<-rownames(rare.otu)
  taxon_table$label<-rownames(taxon_table)
  taxon_newtable<-taxon_table[which(taxon_table$label %in% rare.otu$name),]
  taxon_newtable<-subset(taxon_newtable,select = -label)

  otu_removal<-abun[setdiff(rownames(abun),rownames(rare.otu)),]

  export_object<-list(otu_table=otu_table,tax_table=taxon_newtable,otu_removal=otu_removal)
  cat("\n","totally export",length(export_object),"objects: ")

  cat("\n","otu_table is the abundance table after subsampling.")
  cat("\n","tax_table is the taxonomy table after subsampling.")
  cat("\n","otu_removal is the abundance table of the removal OTUs.","\n")

  return(export_object)

}




"subsampleMicrochat" <- function(tidymchatobj,
                                 filter=TRUE,
                                 sample.size=NULL,
                                 export_path="data") {

  ###subsample samples to uniform size

  if (class(tidymchatobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object
")
  }

  otu_table <- tidymchatobj$otu_table
  taxon_table <- tidymchatobj$taxon_table
  tree <- tidymchatobj$tree
  param_table<-tidymchatobj$param_table


  if (filter) {
    if (length(which(colSums(otu_table)==0)) !=0) otu_table1<-otu_table[,-which(colSums(otu_table)==0)]
    if (length(which(colSums(otu_table)==0)) ==0) otu_table1<-otu_table
    data2<-sub2sample(otudata=otu_table1,taxon_table=taxon_table,sample.size=sample.size)
    otu_tablex<-otu_table
    otu_tablexx<-data2$otu_table[which(rownames(data2$otu_table) %in%
                                         rownames(otu_tablex)) ,]

    newtab<-matrix(0,ncol = ncol(otu_table), nrow = nrow(otu_tablexx))%>%data.frame()
    colnames(newtab)<-colnames(otu_table)
    rownames(newtab)<-rownames(otu_tablexx)

    newtab[,which(colnames(otu_tablex) %in%
                    colnames(data2$otu_table))] <- otu_tablexx

    data2<-list(otu_table=newtab,tax_table=data2$tax_table)
    } else {
    data2<-list(otu_table=otu_table,tax_table=taxon_table)
    }

  otu_table<-data2$otu_table
  taxon_table<-data2$tax_table

  tree_filter<-taxon_filter(tree,taxon_table,otu_table)
  otu_table<-tree_filter$abun
  taxon_table<-tree_filter$taxon
  tree<-tree_filter$tree

  ###data output
  data_export_to_file(otu_table,filename="otu_filter_subsample",export_path)
  data_export_to_file(taxon_table,filename="taxon_filter_subsample",export_path)
  if (!is.null(tree)) write.tree(tree,file = paste(export_path,"/data/tree.nwk",sep = ""))
  if (!is.null(param_table)) data_export_to_file(param_table,filename="param_filter_subsample",export_path)

  ###or use this
  #for (j in names(param_table)) {
  #  newparam<-param_table[[which(names(param_table)==j)]]
  #  data_export_to_file(newparam,filename=j,export_path)
  #}

  ###create S3 object
  submicrochatobj <- list(otu_table,taxon_table,tree,param_table)

  class(submicrochatobj) <- "microchat"

  names(submicrochatobj)[[1]]<-"otu_table"
  names(submicrochatobj)[[2]]<-"taxon_table"
  names(submicrochatobj)[[3]]<-"tree"
  names(submicrochatobj)[[4]]<-"param_table"

  samplenum<-length(unique(substr(colnames(otu_table),start=1,stop=2)))

  cat("\n","otu_table: ",ncol(submicrochatobj[[1]]), " Samples with ",samplenum," groups", sep = "")
  cat("\n","taxon_table: ",nrow(submicrochatobj[[2]]), " OTUs with ", ncol(submicrochatobj[[2]]), " taxa ", sep = "" )
  cat("\n","sample.size selected: ",colSums(otu_table)[1],sep = "","\n")
  print(submicrochatobj[[3]])
  return(submicrochatobj)
}



"data_export_to_file" <- function(otudata,filename="otu_filter_subsample",
                                  filepathX) {
  options(warn = -1)
  if (class(otudata)!="list") {
  otudata$name<-rownames(otudata)
  abun<-otudata
  abun<-cbind(abun$name,abun[,1:(length(colnames(abun))-1)])
  colnames(abun)[1]<-"name"
  filepath0<-paste(filepathX,"/data",sep = "")
  dir.create(filepath0, recursive = TRUE)
  filepath =paste(filepath0,"/",sep = "")
  file1=paste(filepath, filename,".txt",sep = "" )
  write.table(abun,file = file1,  row.names = FALSE,quote = FALSE, sep = "\t")
  } else {
    for (j in names(otudata)) {
      newparam<-otudata[[which(names(otudata)==j)]]
      otudata1<-newparam
      otudata1$name<-rownames(otudata1)
      abun<-otudata1
      abun<-cbind(abun$name,abun[,1:(length(colnames(abun))-1)])
      colnames(abun)[1]<-"name"
      filepath0<-paste(filepathX,"/data",sep = "")
      dir.create(filepath0, recursive = TRUE)
      filepath =paste(filepath0,"/",sep = "")
      file1=paste(filepath, j,".txt",sep = "" )
      write.table(abun,file = file1,  row.names = FALSE,quote = FALSE, sep = "\t")
    }

  }

  cat("\n","-----------------------------------------------------------------------")
  cat("\n","Subsequent analysis used the files located in ",paste(getwd(),"/",filepath,sep = ""))
  cat("\n","-----------------------------------------------------------------------")

}
