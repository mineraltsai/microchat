"rda1Visual"<-function (x,
                        point.shape=point.shape,
                       main.barcolor=c("gray95", "black"),
                       main.pointcolor=c("gray95", "black"),
                       hp.color=c("gray95", "black"),
                       panel.color1=c("white", "gray95"),
                       panel.color2=c("gray95", "black"),
                       plot.hp = FALSE,
                       order.part = "effect",
                       decreasing.part = TRUE,
                       order.var = TRUE,
                       decreasing.var = TRUE,
                       cutoff = -1, nVar = 30,
                       col.width = 0.6,
                       pch.size = 3,
                       line.lwd = 0.5,
                       show.effect = TRUE,
                       effect.cex = 2.7,
                       title.cex = 10,
                       axis.cex = 8,
                       height.ratio = c(2,1),
                       width.ratio = c(1, 3))
{
  Constrained <- 100 * x$Total_explained_variation
  Var.part <- as.data.frame(x$Var.part)
  Var.part <- Var.part[-nrow(Var.part), ]
  Var.part$Var <- gsub("^\\s+|\\s+$", "", gsub("and ", "",
                                               gsub("Common to ", "", gsub("Unique to ", "", rownames(Var.part)))))
  Hier.part <- as.data.frame(x$Hier.part)
  Hier.part$Var <- rownames(Hier.part)
  Var.part$Fractions <- 100 * Var.part$Fractions
  Hier.part$Individual <- 100 * Hier.part$Individual
  Var.part$inter <- lengths(strsplit(Var.part$Var, ", "))
  Var.part$valid <- apply(Var.part[1], 1, function(x) ifelse(x <=
                                                               0, "0", "1"))
  Var.part <- Var.part[which(Var.part$Fractions >= 100 * cutoff),
  ]
  if (order.part == "effect")
    Var.part <- Var.part[order(Var.part$Fractions, decreasing = decreasing.part),
    ] else if (order.part == "degree")
      Var.part <- Var.part[order(Var.part$inter, Var.part$Fractions,
                                 decreasing = c(!decreasing.part, TRUE)), ]



  if (nrow(Var.part) > nVar)
    Var.part <- Var.part[1:nVar, ]
  Var.part$Var <- factor(Var.part$Var, levels = Var.part$Var)

  ysize.max<-max(Var.part$"Fractions")
  ysize.min<-min(Var.part$"Fractions")
  ysize.max.use=(ysize.max%>%as.integer())+1
  if (ysize.min>0) {
    ysize.min.use=0
  } else {
    ysize.min.use=ysize.min
  }

  p.vp <- ggplot2::ggplot(data = Var.part,
                          aes_string(x = "Var", y = "Fractions", fill = "valid")) +
    ggplot2::geom_col(width = col.width) +
    ggplot2::scale_fill_manual(values = main.barcolor,
                               limits = c("0", "1")) +
    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(),
                   axis.line.y = element_line(color = "black"),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(color = "black", size = axis.cex),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_line(color = "black"),
                   axis.title = element_text(color = "black", size = title.cex),
                   plot.title = element_text(hjust = 0.5,size = title.cex),
                   legend.position = "none") +
    ggplot2::labs(y = "Fractions (%)", x = "")+
    ggplot2::coord_cartesian(ylim=c(ysize.min.use,ysize.max.use))



  if (show.effect)
    p.vp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Fractions",
                                                      vjust = ifelse(Var.part$Fractions >= 0, -0.2, 1.2)),
                                           color = "black",
                                           family="serif",
                                           size = effect.cex)
  p.vp <- p.vp + ggplot2::geom_hline(yintercept = 0)




  Fractions <- NULL
  for (i in rownames(Hier.part)) Fractions <- c(Fractions,
                                                sum(Var.part[grep(i, Var.part$Var), "Fractions"]))
  Var.exp <- data.frame(Var = rownames(Hier.part), Fractions = Fractions)
  Var.exp$valid <- apply(Var.exp[2], 1, function(x) ifelse(x <=
                                                             0, "0", "1"))
  Var.exp <- Var.exp[which(Var.exp$Fractions >= 100 * cutoff),
  ]
  if (order.var)
    Var.exp <- Var.exp[order(Var.exp$Fractions, decreasing = !decreasing.var),
    ]
  Var.exp$Var <- factor(Var.exp$Var, levels = Var.exp$Var)

  p.exp <- ggplot2::ggplot(data = Var.exp, aes_string(x = "Var",
                                                      y = "Fractions", fill = "valid")) + ggplot2::geom_col(width = col.width) +
    ggplot2::scale_fill_manual(values = main.pointcolor,
                               limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                      panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                      axis.text.x = element_text(color = "black", size = axis.cex),
                                                                      axis.text.y = element_blank(), axis.ticks.x = element_line(color = "black"),
                                                                      axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                size = title.cex), legend.position = "none") + ggplot2::coord_flip() +
    ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                         ifelse(min(Var.exp$Fractions) < 0, 0.3, 0)))) +
    ggplot2::labs(y = "Fractions (%)", x = NULL)

  if (plot.hp) {
    Hier.part <- Hier.part[which(Hier.part$Individual >=
                                   100 * cutoff), ]
    if (order.var)
      Hier.part <- Hier.part[order(Hier.part$Individual,
                                   decreasing = !decreasing.var), ]
    Hier.part$Var <- factor(Hier.part$Var, levels = Hier.part$Var)
    Hier.part$valid <- apply(Hier.part[3], 1, function(x) ifelse(x <=
                                                                   0, "0", "1"))
    p.hp <- ggplot2::ggplot(data = Hier.part, aes_string(x = "Var",
                                                         y = "Individual", fill = "valid")) + ggplot2::geom_col(width = col.width) +
      ggplot2::scale_fill_manual(values = hp.color, limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                                             panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                                             axis.text.x = element_text(color = "black", size = axis.cex),
                                                                                             axis.text.y = element_blank(), axis.ticks = element_line(color = "black"),
                                                                                             axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                                       size = title.cex), legend.position = "none") +
      ggplot2::coord_flip() + ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                                                   ifelse(min(Hier.part$Individual) < 0, 0.3, 0)))) +
      ggplot2::labs(y = "Individual (%)", x = NULL)
    if (show.effect)
      p.hp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Individual",
                                                        hjust = ifelse(Hier.part$Individual >= 0, 1.2,
                                                                       -0.2)), color = "black", size = effect.cex)
    p.hp <- p.hp + ggplot2::geom_hline(yintercept = 0)
  }
  panel <- NULL
  if (plot.hp)
    Var <- Hier.part else Var <- Var.exp
  for (i in Var$Var) for (j in Var.part$Var) panel <- rbind(panel,
                                                            c(i, j, 0))
  panel <- data.frame(panel, stringsAsFactors = FALSE)
  panel <- panel[which(panel$X2 %in% Var.part$Var), ]
  panel$X1 <- factor(panel$X1, levels = Var$Var)
  panel$X2 <- factor(panel$X2, levels = Var.part$Var)
  for (i in 1:nrow(Var.part)) {
    i <- as.character(Var.part[i, "Var"])
    for (j in unlist(strsplit(i, ", "))) panel[which(panel$X1 ==
                                                       j & panel$X2 == i), "X3"] <- "1"
  }
  panel$X4 <- apply(panel, 1, function(x) ifelse(which(levels(panel$X1) ==
                                                         x[1])%%2 == 0, "0", "1"))


  Var.exp1<-left_join(panel,Hier.part,by=c("X1"="Var"))
  Var.exp1$yaxis<-Var.exp1$X1
  Var.exp1<- Var.exp1[order(Var.exp1$Individual),]
  Var.exp1$X1 <- factor(Var.exp1$X1, levels = unique(Var.exp1$X1))

  Var.exp1$X4 <- apply(Var.exp1, 1, function(x) ifelse(which(levels(Var.exp1$X1) ==
                                                               x[1])%%2 == 0, "0", "1"))

  Hier.part1<-Hier.part
  Hier.part1<-Hier.part1[order(Hier.part1$Individual,decreasing = FALSE),]
  Hier.part1$Var<-factor(Hier.part1$Individual,levels = unique(Hier.part1$Individual))

  p3 <- ggplot(data = Hier.part1,
               aes(x = Var, y = Individual)) +
    geom_text(aes(y = 0,x = Var,
                  label = scales::percent(Individual,
                                          scale = 1,
                                          suffix="",
                                          accuracy = 0.01)),
              family = "serif",
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          text = element_text(family = "serif"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())


  if (length(point.shape)==1) {
    point.shape<-c(point.shape,point.shape)
  } else if (length(point.shape)==2) {
    point.shape<-point.shape
  } else if (length(point.shape)>2) {
    stop("Please provide no more than two shapes !!!")
  }
  point.shape.use<-point.shape
  names(point.shape.use)<-unique(Var.exp1$X3)[order(unique(Var.exp1$X3))]


  p.panel <- ggplot2::ggplot(data = Var.exp1,
                             aes(x = X2, y = yaxis, color = X3, fill = X4,shape=X3)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::geom_point(size = pch.size) +
    ggplot2::scale_fill_manual(values = panel.color1,
                               limits = c("0", "1")) +
    ggplot2::scale_shape_manual(values = point.shape.use)+
    ggplot2::scale_color_manual(values = panel.color2,
                                limits = c("0", "1")) +
    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(),
                   axis.text = element_text(color = "black",
                                            size = title.cex),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::labs(y = NULL, x = NULL)

  for (i in levels(Var.exp1$X2)) {
    panel.i <- Var.exp1[which(Var.exp1$X2 == i & Var.exp1$X3 == 1),]
    panel.i <- panel.i[order(panel.i$Individual), ]
    if (nrow(panel.i) > 1)
    {
      p.panel <- p.panel +
        ggplot2::annotate("segment",
                          x = i, xend = i,
                          y = panel.i[1, 1], yend = panel.i[nrow(panel.i),1],
                          color = "black", size = line.lwd)
    }
  }

  p.vp<-p.vp+theme(text=element_text(family = "serif"))
  p.panel<-p.panel+theme(text=element_text(family = "serif"))


  if (plot.hp)
  {
    p.vp<-p.vp +
    annotate("text",  x = 19, y = max(Var.part$Fractions)*1.05, family="serif",
             label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))

  addp<-ggplot()+theme(panel.background = element_rect(fill = "white"))
  xxp<-p.panel+p3+ patchwork::plot_layout(widths  = c(9.5,0.5),byrow=TRUE)
  xxpx<-p.vp+addp+ patchwork::plot_layout(widths  = c(9.5,0.5),byrow=TRUE)

  (xxpx + p.panel  + patchwork::plot_layout(heights = c(6.5,3.5),byrow=TRUE))+
    p3

  } else {
    p.vp<-p.vp+ scale_y_reverse()+
      annotate("text",  x = 19, y = max(Var.part$Fractions)*1.05, family="serif",
               label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))
    xxp<-p.panel+p3+ patchwork::plot_layout(widths = c(9.5,0.5),byrow=TRUE)
    xxp+ p.vp   + patchwork::plot_layout(heights = c(3.5,6.5),byrow=TRUE)
  }


}

"rda2Visual"<-function (x,
                        point.shape=point.shape,
                        main.barcolor=c("gray95", "black"),
                        main.pointcolor=c("gray95", "black"),
                        hp.color=c("gray95", "black"),
                        panel.color1=c("white", "gray95"),
                        panel.color2=c("gray95", "black"),
                        plot.hp = FALSE,
                        order.part = "effect",
                        decreasing.part = TRUE,
                        order.var = TRUE,
                        decreasing.var = TRUE,
                        cutoff = -1, nVar = 30,
                        col.width = 0.6,
                        pch.size = 3,
                        line.lwd = 0.5,
                        show.effect = TRUE,
                        effect.cex = 2.7,
                        title.cex = 10,
                        axis.cex = 8,
                        height.ratio = c(2,1),
                        width.ratio = c(1, 3))
{

  Constrained <- 100 * x$Total_explained_variation
  Var.part <- as.data.frame(x$Var.part)
  Var.part <- Var.part[-nrow(Var.part), ]
  Var.part$Var <- gsub("^\\s+|\\s+$", "", gsub("and ", "",
                                               gsub("Common to ", "", gsub("Unique to ", "", rownames(Var.part)))))
  Hier.part <- as.data.frame(x$Hier.part)
  Hier.part$Var <- rownames(Hier.part)
  Var.part$Fractions <- 100 * Var.part$Fractions
  Hier.part$Individual <- 100 * Hier.part$Individual
  Var.part$inter <- lengths(strsplit(Var.part$Var, ", "))
  Var.part$valid <- apply(Var.part[1], 1, function(x) ifelse(x <=
                                                               0, "0", "1"))
  Var.part <- Var.part[which(Var.part$Fractions >= 100 * cutoff),]

  if (order.part == "effect")
  {
    Var.part <- Var.part[order(Var.part$Fractions, decreasing = decreasing.part),]
  } else if (order.part == "degree") {
    Var.part <- Var.part[order(Var.part$inter, Var.part$Fractions,
                               decreasing = c(!decreasing.part, TRUE)), ]
  }
  if (nrow(Var.part) > nVar) {
    Var.part <- Var.part[1:nVar, ]
  }
  Var.part$Var <- factor(Var.part$Var, levels = Var.part$Var)

  ysize.max<-max(Var.part$"Fractions")
  ysize.min<-min(Var.part$"Fractions")
  ysize.max.use=(ysize.max%>%as.integer())+1
  if (ysize.min>0) {
    ysize.min.use=0
  } else {
    ysize.min.use=ysize.min
  }

  p.vp <- ggplot2::ggplot(data = Var.part,
                          aes_string(x = "Var", y = "Fractions", fill = "valid")) +
    ggplot2::geom_col(width = col.width) +
    ggplot2::scale_fill_manual(values = main.barcolor,
                               limits = c("0", "1")) +
    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(),
                   axis.line.y = element_line(color = "black"),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(color = "black", size = axis.cex),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_line(color = "black"),
                   axis.title = element_text(color = "black", size = title.cex),
                   plot.title = element_text(hjust = 0.5,size = title.cex),
                   legend.position = "none") +
    ggplot2::labs(y = "Fractions (%)", x = "")+
    ggplot2::coord_cartesian(ylim=c(ysize.max.use,ysize.min.use))


  if (show.effect)
    p.vp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Fractions",
                                                      vjust = ifelse(Var.part$Fractions >= 0, 1, -0.2)),
                                           color = "black",
                                           family="serif",
                                           size = effect.cex)

  p.vp <- p.vp + ggplot2::geom_hline(yintercept = 0)
  Fractions <- NULL
  for (i in rownames(Hier.part)) {
    Fractions <- c(Fractions,sum(Var.part[grep(i, Var.part$Var), "Fractions"]))
  }

  Var.exp <- data.frame(Var = rownames(Hier.part), Fractions = Fractions)
  Var.exp$valid <- apply(Var.exp[2], 1, function(x) ifelse(x <= 0, "0", "1"))
  Var.exp <- Var.exp[which(Var.exp$Fractions >= 100 * cutoff),]

  if (order.var)
  {
    Var.exp <- Var.exp[order(Var.exp$Fractions, decreasing = !decreasing.var),]
  }

  Var.exp$Var <- factor(Var.exp$Var, levels = Var.exp$Var)

  p.exp <- ggplot2::ggplot(data = Var.exp, aes_string(x = "Var",
                                                      y = "Fractions", fill = "valid")) + ggplot2::geom_col(width = col.width) +
    ggplot2::scale_fill_manual(values = main.pointcolor,
                               limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                      panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                      axis.text.x = element_text(color = "black", size = axis.cex),
                                                                      axis.text.y = element_blank(), axis.ticks.x = element_line(color = "black"),
                                                                      axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                size = title.cex), legend.position = "none") + ggplot2::coord_flip() +
    ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                         ifelse(min(Var.exp$Fractions) < 0, 0.3, 0)))) +
    ggplot2::labs(y = "Fractions (%)", x = NULL)



  if (plot.hp) {
    Hier.part <- Hier.part[which(Hier.part$Individual >=
                                   100 * cutoff), ]
    if (order.var)
      Hier.part <- Hier.part[order(Hier.part$Individual,
                                   decreasing = !decreasing.var), ]
    Hier.part$Var <- factor(Hier.part$Var, levels = Hier.part$Var)
    Hier.part$valid <- apply(Hier.part[3], 1, function(x) ifelse(x <=
                                                                   0, "0", "1"))
    p.hp <- ggplot2::ggplot(data = Hier.part, aes_string(x = "Var",
                                                         y = "Individual", fill = "valid")) + ggplot2::geom_col(width = col.width) +
      ggplot2::scale_fill_manual(values = hp.color, limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                                           panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                                           axis.text.x = element_text(color = "black", size = axis.cex),
                                                                                           axis.text.y = element_blank(), axis.ticks = element_line(color = "black"),
                                                                                           axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                                     size = title.cex), legend.position = "none") +
      ggplot2::coord_flip() + ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                                                   ifelse(min(Hier.part$Individual) < 0, 0.3, 0)))) +
      ggplot2::labs(y = "Individual (%)", x = NULL)
    if (show.effect)
      p.hp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Individual",
                                                        hjust = ifelse(Hier.part$Individual >= 0, 1.2,
                                                                       -0.2)), color = "black", size = effect.cex)
    p.hp <- p.hp + ggplot2::geom_hline(yintercept = 0)
  }
  panel <- NULL
  if (plot.hp) {Var <- Hier.part} else {Var <- Var.exp}
  for (i in Var$Var) for (j in Var.part$Var) panel <- rbind(panel,
                                                            c(i, j, 0))
  panel <- data.frame(panel, stringsAsFactors = FALSE)
  panel <- panel[which(panel$X2 %in% Var.part$Var), ]
  panel$X1 <- factor(panel$X1, levels = Var$Var)
  panel$X2 <- factor(panel$X2, levels = Var.part$Var)
  for (i in 1:nrow(Var.part)) {
    i <- as.character(Var.part[i, "Var"])
    for (j in unlist(strsplit(i, ", "))) panel[which(panel$X1 ==
                                                       j & panel$X2 == i), "X3"] <- "1"
  }
  panel$X4 <- apply(panel, 1, function(x) ifelse(which(levels(panel$X1) ==
                                                         x[1])%%2 == 0, "0", "1"))


  Var.exp1<-left_join(panel,Hier.part,by=c("X1"="Var"))
  Var.exp1$yaxis<-Var.exp1$X1
  Var.exp1<- Var.exp1[order(Var.exp1$Individual,decreasing = TRUE),]
  Var.exp1$X1 <- factor(Var.exp1$X1, levels = unique(Var.exp1$X1))

  Var.exp1$X4 <- apply(Var.exp1, 1, function(x) ifelse(which(levels(Var.exp1$X1) ==
                                                               x[1])%%2 == 0, "0", "1"))

  #Var.exp1$X1<-factor(Var.exp1$X1, levels = unique(Var.exp1$X1)[length(unique(Var.exp1$X1)):1] )

  if (length(point.shape)==1) {
    point.shape<-c(point.shape,point.shape)
  } else if (length(point.shape)==2) {
    point.shape<-point.shape
  } else if (length(point.shape)>2) {
    stop("Please provide no more than two shapes !!!")
  }

  point.shape.use<-point.shape
  names(point.shape.use)<-unique(Var.exp1$X3)[order(unique(Var.exp1$X3))]

  p.panel <- ggplot2::ggplot(data = Var.exp1,
                             aes(x = X2, y = X1, color = X3, fill = X4,shape=X3)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::geom_point(size = pch.size) +
    ggplot2::scale_shape_manual(values = point.shape.use)+
    ggplot2::scale_fill_manual(values = panel.color1,
                               limits = c("0", "1")) +
    ggplot2::scale_color_manual(values = panel.color2,
                                limits = c("0", "1")) +
    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(),
                   axis.text = element_text(color = "black",
                                            size = title.cex),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::labs(y = NULL, x = NULL)

  Hier.part1<-Hier.part
  Hier.part1<-Hier.part1[order(Hier.part1$Individual,decreasing = TRUE),]
  Hier.part1$Var<-factor(Hier.part1$Individual,levels = unique(Hier.part1$Individual))

  p3 <- ggplot(data = Hier.part1,
               aes(x = Var, y = Individual)) +
    geom_text(aes(y = 0,x = Var,
                  label = scales::percent(Individual,
                                          scale = 1,
                                          suffix="",
                                          accuracy = 0.01)),
              family = "serif",
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          text = element_text(family = "serif"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())



  for (i in levels(Var.exp1$X2)) {
    panel.i <- Var.exp1[which(Var.exp1$X2 == i & Var.exp1$X3 == 1),
    ]
    panel.i <- panel.i[order(panel.i$X1), ]
    if (nrow(panel.i) > 1)
      p.panel <- p.panel + ggplot2::annotate("segment",
                                             x = i, xend = i, y = panel.i[1, 1], yend = panel.i[nrow(panel.i),
                                                                                                1], color = "black", size = line.lwd)
  }


  p.vp<-p.vp+theme(text=element_text(family = "serif"))
  p.panel<-p.panel+theme(text=element_text(family = "serif"))


  if (plot.hp)
  {
    p.vp<-p.vp +
    annotate("text",  x = 19, y = max(Var.part$Fractions)*1.05, family="serif",
             label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))

    xp<-p.panel +p3 +patchwork::plot_layout(widths=c(9.5,0.5),byrow=TRUE)
    p.vp +xp+  + patchwork::plot_layout(heights = c(3.5,6.5),byrow=TRUE)
  } else {
    p.vp<-p.vp+ scale_y_reverse()+
      annotate("text",  x = 19, y = max(Var.part$Fractions)*1.05, family="serif",
               label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))

    xp<-p.panel +p3 +patchwork::plot_layout(widths=c(9.5,0.5),byrow=TRUE)
    xp+ p.vp   + patchwork::plot_layout(heights = c(3.5,6.5),byrow=TRUE)
  }


}


"rda3Visual"<-function (x,
                        point.shape=point.shape,
                        main.barcolor=c("gray95", "black"),
                        main.pointcolor=c("gray95", "black"),
                        hp.color=c("gray95", "black"),
                        panel.color1=c("white", "gray95"),
                        panel.color2=c("gray95", "black"),
                        plot.hp = FALSE,
                        order.part = "effect",
                        decreasing.part = TRUE,
                        order.var = TRUE,
                        decreasing.var = TRUE,
                        cutoff = -1, nVar = 30,
                        col.width = 0.6,
                        pch.size = 3,
                        line.lwd = 0.5,
                        show.effect = TRUE,
                        effect.cex = 2.7,
                        title.cex = 10,
                        axis.cex = 8,
                        height.ratio = c(2,1),
                        width.ratio = c(1, 3))
{
  Constrained <- 100 * x$Total_explained_variation
  Var.part <- as.data.frame(x$Var.part)
  Var.part <- Var.part[-nrow(Var.part), ]
  Var.part$Var <- gsub("^\\s+|\\s+$", "", gsub("and ", "",
                                               gsub("Common to ", "", gsub("Unique to ", "", rownames(Var.part)))))
  Hier.part <- as.data.frame(x$Hier.part)
  Hier.part$Var <- rownames(Hier.part)
  Var.part$Fractions <- 100 * Var.part$Fractions
  Hier.part$Individual <- 100 * Hier.part$Individual
  Var.part$inter <- lengths(strsplit(Var.part$Var, ", "))
  Var.part$valid <- apply(Var.part[1], 1, function(x) ifelse(x <=
                                                               0, "0", "1"))
  Var.part <- Var.part[which(Var.part$Fractions >= 100 * cutoff),
  ]
  if (order.part == "effect")
    Var.part <- Var.part[order(Var.part$Fractions, decreasing = decreasing.part),
    ] else if (order.part == "degree")
      Var.part <- Var.part[order(Var.part$inter, Var.part$Fractions,
                                 decreasing = c(!decreasing.part, TRUE)), ]



  if (nrow(Var.part) > nVar)
    Var.part <- Var.part[1:nVar, ]
  Var.part$Var <- factor(Var.part$Var, levels = Var.part$Var[length(Var.part$Var):1])

  ysize.max<-max(Var.part$"Fractions")
  ysize.min<-min(Var.part$"Fractions")
  ysize.max.use=(ysize.max%>%as.integer())+1
  if (ysize.min>0) {
    ysize.min.use=0
  } else {
    ysize.min.use=ysize.min
  }

  p.vp <-
    ggplot2::ggplot(data = Var.part,
                    aes_string(x = "Var", y = "Fractions",
                               fill = "valid")) +
    ggplot2::geom_col(width = col.width) +
    ggplot2::scale_fill_manual(values = main.barcolor,
                               limits = c("0", "1")) +

    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(),
                   axis.line.y = element_line(color = "black"),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(color = "black",
                                              size = axis.cex),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_line(color = "black"),
                   axis.title = element_text(color = "black", size = title.cex),
                   plot.title = element_text(hjust = 0.5,size = title.cex),
                   legend.position = "none") +
    ggplot2::labs(y = "Fractions (%)", x = "")+
    ggplot2::theme(axis.line.y.left  = element_line(size = NA),
                   axis.title.x = element_blank(),
                   axis.ticks.y.left = element_blank(),
                   axis.title.y.left= element_blank(),
                   axis.text.y.left = element_blank())+
    ggplot2::scale_y_continuous(sec.axis = sec_axis(~.*1))+
    ggplot2::coord_cartesian(ylim=c(ysize.min.use,ysize.max.use))

  if (show.effect) p.vp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Fractions",
                                                      vjust = ifelse(Var.part$Fractions >= 0, -0.2, 1.2)),
                                           color = "black",
                                           family="serif",
                                           size = effect.cex)

  p.vp <- p.vp + ggplot2::geom_hline(yintercept = 0)
  Fractions <- NULL
  for (i in rownames(Hier.part)) Fractions <- c(Fractions,
                                                sum(Var.part[grep(i, Var.part$Var), "Fractions"]))
  Var.exp <- data.frame(Var = rownames(Hier.part), Fractions = Fractions)
  Var.exp$valid <- apply(Var.exp[2], 1, function(x) ifelse(x <=
                                                             0, "0", "1"))
  Var.exp <- Var.exp[which(Var.exp$Fractions >= 100 * cutoff),]
  if (order.var)
    Var.exp <- Var.exp[order(Var.exp$Fractions, decreasing = !decreasing.var),]
  Var.exp$Var <- factor(Var.exp$Var, levels = Var.exp$Var)

  p.exp <- ggplot2::ggplot(data = Var.exp, aes_string(x = "Var",
                                                      y = "Fractions", fill = "valid")) + ggplot2::geom_col(width = col.width) +
    ggplot2::scale_fill_manual(values = main.pointcolor,
                               limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                      panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                      axis.text.x = element_text(color = "black", size = axis.cex),
                                                                      axis.text.y = element_blank(), axis.ticks.x = element_line(color = "black"),
                                                                      axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                size = title.cex), legend.position = "none") + ggplot2::coord_flip() +
    ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                         ifelse(min(Var.exp$Fractions) < 0, 0.3, 0)))) +
    ggplot2::labs(y = "Fractions (%)", x = NULL)



  if (plot.hp) {
    Hier.part <- Hier.part[which(Hier.part$Individual >=
                                   100 * cutoff), ]
    if (order.var)
      Hier.part <- Hier.part[order(Hier.part$Individual,
                                   decreasing = !decreasing.var), ]
    Hier.part$Var <- factor(Hier.part$Var, levels = Hier.part$Var)
    Hier.part$valid <- apply(Hier.part[3], 1, function(x) ifelse(x <=
                                                                   0, "0", "1"))
    p.hp <- ggplot2::ggplot(data = Hier.part, aes_string(x = "Var",
                                                         y = "Individual", fill = "valid")) + ggplot2::geom_col(width = col.width) +
      ggplot2::scale_fill_manual(values = hp.color, limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                                           panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                                           axis.text.x = element_text(color = "black", size = axis.cex),
                                                                                           axis.text.y = element_blank(), axis.ticks = element_line(color = "black"),
                                                                                           axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                                     size = title.cex), legend.position = "none") +
      ggplot2::coord_flip() + ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                                                   ifelse(min(Hier.part$Individual) < 0, 0.3, 0)))) +
      ggplot2::labs(y = "Individual (%)", x = NULL)
   if (show.effect)
      p.hp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Individual",
                                                        hjust = ifelse(Hier.part$Individual >= 0, 1.2,
                                                                       -0.2)), color = "black", size = effect.cex)
    p.hp <- p.hp + ggplot2::geom_hline(yintercept = 0)
  }
  panel <- NULL
  if (plot.hp)
    Var <- Hier.part else Var <- Var.exp
  for (i in Var$Var) for (j in Var.part$Var) panel <- rbind(panel,
                                                            c(i, j, 0))
  panel <- data.frame(panel, stringsAsFactors = FALSE)
  panel <- panel[which(panel$X2 %in% Var.part$Var), ]
  panel$X1 <- factor(panel$X1, levels = Var$Var)
  panel$X2 <- factor(panel$X2, levels = Var.part$Var)
  for (i in 1:nrow(Var.part)) {
    i <- as.character(Var.part[i, "Var"])
    for (j in unlist(strsplit(i, ", "))) panel[which(panel$X1 ==
                                                       j & panel$X2 == i), "X3"] <- "1"
  }
  panel$X4 <- apply(panel, 1, function(x) ifelse(which(levels(panel$X1) ==
                                                         x[1])%%2 == 0, "0", "1"))


  Var.exp1<-left_join(panel,Hier.part,by=c("X1"="Var"))
  Var.exp1$yaxis<-Var.exp1$X1
  Var.exp1<- Var.exp1[order(Var.exp1$Individual),]
  Var.exp1$X1 <- factor(Var.exp1$X1, levels = unique(Var.exp1$X1))

  Var.exp1$X4 <- apply(Var.exp1, 1, function(x) ifelse(which(levels(Var.exp1$X1) ==
                                                               x[1])%%2 == 0, "0", "1"))

  Hier.part1<-Hier.part
  Hier.part1<-Hier.part1[order(Hier.part1$Individual,decreasing = FALSE),]
  Hier.part1$Var<-factor(Hier.part1$Individual,levels = unique(Hier.part1$Individual))

  p3 <- ggplot(data = Hier.part1,
               aes(x = Var, y = Individual)) +
    geom_text(aes(y = 0,x = Var,
                  label = scales::percent(Individual,
                                          scale = 1,
                                          suffix="",
                                          accuracy = 0.01)),
              family = "serif",
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          text = element_text(family = "serif"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())

  Hier.part2<-tibble::rownames_to_column(Hier.part1,var = "X1")
  Hier.part2$X1<-factor(Hier.part2$X1,levels = unique(Hier.part2$X1)[length(Hier.part2$X1):1])
  p4 <- ggplot(data = Hier.part2,
               aes(x = Var, y = Individual)) +
    geom_text(aes(y = 0,x = Var,
                  label = X1),color="white",fill="white",
              family = "serif",
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          text = element_text(family = "serif"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())

  Var.exp1$X2 <- factor(Var.exp1$X2, levels = unique(Var.exp1$X2[length(Var.exp1$X2):1]))
  Var.exp1<-Var.exp1[order(Var.exp1$Individual,decreasing = TRUE),]
  unique(Var.exp1$X1[length(Var.exp1$X1):1])
  Var.exp1$X1 <- factor(Var.exp1$X1, levels = unique(Var.exp1$X1[length(Var.exp1$X1):1]))


  if (length(point.shape)==1) {
    point.shape<-c(point.shape,point.shape)
  } else if (length(point.shape)==2) {
    point.shape<-point.shape
  } else if (length(point.shape)>2) {
    stop("Please provide no more than two shapes !!!")
  }

  point.shape.use<-point.shape
  names(point.shape.use)<-unique(Var.exp1$X3)[order(unique(Var.exp1$X3))]

  p.panel <- ggplot2::ggplot(data = Var.exp1,
                             aes(x = X2, y = X1, color = X3, fill = X4,shape=X3)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::geom_point(size = pch.size) +
    ggplot2::scale_shape_manual(values = point.shape.use)+
    ggplot2::scale_fill_manual(values = panel.color1,
                               limits = c("0", "1")) +
    ggplot2::scale_color_manual(values = panel.color2,
                                limits = c("0", "1")) +
    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(),
                   axis.text = element_text(color = "black",
                                            size = title.cex),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::labs(y = NULL, x = NULL)+
    ggplot2::theme(axis.text.y.left = element_blank())




  for (i in levels(Var.exp1$X2)) {
    panel.i <- Var.exp1[which(Var.exp1$X2 == i & Var.exp1$X3 == 1),]
    panel.i <- panel.i[order(panel.i$Individual), ]
    if (nrow(panel.i) > 1)
    {
      p.panel <- p.panel +
        ggplot2::annotate("segment",
                          x = i, xend = i,
                          y = panel.i[1, 1], yend = panel.i[nrow(panel.i),1],
                          color = "black", size = line.lwd)
    }
  }

  p.vp<-p.vp+theme(text=element_text(family = "serif"))
  p.panel<-p.panel+theme(text=element_text(family = "serif"))


  if (plot.hp)
  {
    p.vp<-p.vp +
      annotate("text",  x = 11, y = max(Var.part$Fractions)*1.05, family="serif",
               label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))

    addp<-ggplot(data = Hier.part2,
                 aes(x = Var, y = Individual)) +
      geom_text(aes(y = 0,x = Var,
                    label = X1),color="white",
                family = "serif",
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_blank(),axis.ticks = element_blank())


      (addp+p.vp+patchwork::plot_layout(widths = c(0.5,9.5)))/
      (p3+p.panel+patchwork::plot_layout(widths = c(0.5,9.5)))+
      patchwork::plot_layout(heights = c(6.5,3.5))


  } else {
    p.vp<-p.vp+ scale_y_reverse()+
      annotate("text",  x = 25, y = max(Var.part$Fractions)*1.05, family="serif",
               label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))
    xxp<-p.panel+p3+ patchwork::plot_layout(widths = c(9.5,0.5),byrow=TRUE)
    xxp+ p.vp   + patchwork::plot_layout(heights = c(3.5,6.5),byrow=TRUE)
  }
}




"rda4Visual"<-function (x,
                        point.shape=point.shape,
                        main.barcolor=c("gray95", "black"),
                        main.pointcolor=c("gray95", "black"),
                        hp.color=c("gray95", "black"),
                        panel.color1=c("white", "gray95"),
                        panel.color2=c("gray95", "black"),
                        plot.hp = FALSE,
                        order.part = "effect",
                        decreasing.part = TRUE,
                        order.var = TRUE,
                        decreasing.var = TRUE,
                        cutoff = -1, nVar = 30,
                        col.width = 0.6,
                        pch.size = 3,
                        line.lwd = 0.5,
                        show.effect = TRUE,
                        effect.cex = 2.7,
                        title.cex = 10,
                        axis.cex = 8,
                        height.ratio = c(2,1),
                        width.ratio = c(1, 3))
{

  Constrained <- 100 * x$Total_explained_variation
  Var.part <- as.data.frame(x$Var.part)
  Var.part <- Var.part[-nrow(Var.part), ]
  Var.part$Var <- gsub("^\\s+|\\s+$", "", gsub("and ", "",
                                               gsub("Common to ", "", gsub("Unique to ", "", rownames(Var.part)))))
  Hier.part <- as.data.frame(x$Hier.part)
  Hier.part$Var <- rownames(Hier.part)
  Var.part$Fractions <- 100 * Var.part$Fractions
  Hier.part$Individual <- 100 * Hier.part$Individual
  Var.part$inter <- lengths(strsplit(Var.part$Var, ", "))
  Var.part$valid <- apply(Var.part[1], 1, function(x) ifelse(x <=
                                                               0, "0", "1"))
  Var.part <- Var.part[which(Var.part$Fractions >= 100 * cutoff),]

  if (order.part == "effect") {
    Var.part <- Var.part[order(Var.part$Fractions, decreasing = decreasing.part),]
  } else if (order.part == "degree") {
    Var.part <- Var.part[order(Var.part$inter, Var.part$Fractions,
                               decreasing = c(!decreasing.part, TRUE)), ]
  }
  if (nrow(Var.part) > nVar) Var.part <- Var.part[1:nVar, ]
    Var.part$Var <- factor(Var.part$Var, levels = Var.part$Var[length(Var.part$Var):1])

    ysize.max<-max(Var.part$"Fractions")
    ysize.min<-min(Var.part$"Fractions")
    ysize.max.use=(ysize.max%>%as.integer())+1
    if (ysize.min>0) {
      ysize.min.use=0
    } else {
      ysize.min.use=ysize.min
    }

    p.vp <-
      ggplot2::ggplot(data = Var.part,
                      aes_string(x = "Var", y = "Fractions",
                                 fill = "valid")) +
      ggplot2::geom_col(width = col.width) +
      ggplot2::scale_fill_manual(values = main.barcolor,
                                 limits = c("0", "1")) +

      ggplot2::theme(panel.grid = element_blank(),
                     panel.background = element_blank(),
                     axis.line.y = element_line(color = "black"),
                     axis.text.x = element_blank(),
                     axis.text.y = element_text(color = "black",
                                                size = axis.cex),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_line(color = "black"),
                     axis.title = element_text(color = "black", size = title.cex),
                     plot.title = element_text(hjust = 0.5,size = title.cex),
                     legend.position = "none") +
      ggplot2::labs(y = "Fractions (%)", x = "")+
      ggplot2::theme(axis.line.y.left  = element_line(size = NA),
                     axis.title.x = element_blank(),
                     axis.ticks.y.left = element_blank(),
                     axis.title.y.left= element_blank(),
                     axis.text.y.left = element_blank())+
      ggplot2::scale_y_continuous(sec.axis = sec_axis(~.*1))+
      ggplot2::coord_cartesian(ylim=c(ysize.max.use,ysize.min.use))

  if (show.effect)
    p.vp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Fractions",
                                                      vjust = ifelse(Var.part$Fractions >= 0, 1, 1.2)),
                                           color = "black",
                                           family="serif",
                                           size = effect.cex)

  p.vp <- p.vp + ggplot2::geom_hline(yintercept = 0)
  Fractions <- NULL
  for (i in rownames(Hier.part)) {
    Fractions <- c(Fractions,sum(Var.part[grep(i, Var.part$Var), "Fractions"]))
  }

  Var.exp <- data.frame(Var = rownames(Hier.part), Fractions = Fractions)
  Var.exp$valid <- apply(Var.exp[2], 1, function(x) ifelse(x <= 0, "0", "1"))
  Var.exp <- Var.exp[which(Var.exp$Fractions >= 100 * cutoff),]

  if (order.var) {
    Var.exp <- Var.exp[order(Var.exp$Fractions, decreasing = !decreasing.var),]
  }

  Var.exp$Var <- factor(Var.exp$Var, levels = Var.exp$Var)

  p.exp <- ggplot2::ggplot(data = Var.exp, aes_string(x = "Var",
                                                      y = "Fractions", fill = "valid")) + ggplot2::geom_col(width = col.width) +
    ggplot2::scale_fill_manual(values = main.pointcolor,
                               limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                      panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                      axis.text.x = element_text(color = "black", size = axis.cex),
                                                                      axis.text.y = element_blank(), axis.ticks.x = element_line(color = "black"),
                                                                      axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                size = title.cex), legend.position = "none") + ggplot2::coord_flip() +
    ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                         ifelse(min(Var.exp$Fractions) < 0, 0.3, 0)))) +
    ggplot2::labs(y = "Fractions (%)", x = NULL)



  if (plot.hp) {
    Hier.part <- Hier.part[which(Hier.part$Individual >=
                                   100 * cutoff), ]
    if (order.var)
      Hier.part <- Hier.part[order(Hier.part$Individual,
                                   decreasing = !decreasing.var), ]
    Hier.part$Var <- factor(Hier.part$Var, levels = Hier.part$Var)
    Hier.part$valid <- apply(Hier.part[3], 1, function(x) ifelse(x <=
                                                                   0, "0", "1"))
    p.hp <- ggplot2::ggplot(data = Hier.part, aes_string(x = "Var",
                                                         y = "Individual", fill = "valid")) + ggplot2::geom_col(width = col.width) +
      ggplot2::scale_fill_manual(values = hp.color, limits = c("0", "1")) + ggplot2::theme(panel.grid = element_blank(),
                                                                                           panel.background = element_blank(), axis.line.x = element_line(color = "black"),
                                                                                           axis.text.x = element_text(color = "black", size = axis.cex),
                                                                                           axis.text.y = element_blank(), axis.ticks = element_line(color = "black"),
                                                                                           axis.ticks.y = element_blank(), axis.title = element_text(color = "black",
                                                                                                                                                     size = title.cex), legend.position = "none") +
      ggplot2::coord_flip() + ggplot2::scale_y_reverse(expand = expansion(mult = c(0.3,
                                                                                   ifelse(min(Hier.part$Individual) < 0, 0.3, 0)))) +
      ggplot2::labs(y = "Individual (%)", x = NULL)
    if (show.effect)
      p.hp$layers[[2]] <- ggplot2::geom_text(aes_string(label = "Individual",
                                                        hjust = ifelse(Hier.part$Individual >= 0, 1.2,
                                                                       -0.2)), color = "black", size = effect.cex)

    p.hp <- p.hp + ggplot2::geom_hline(yintercept = 0)
  }
  panel <- NULL
  if (plot.hp) {Var <- Hier.part} else {Var <- Var.exp}
  for (i in Var$Var) for (j in Var.part$Var) panel <- rbind(panel,
                                                            c(i, j, 0))
  panel <- data.frame(panel, stringsAsFactors = FALSE)
  panel <- panel[which(panel$X2 %in% Var.part$Var), ]
  panel$X1 <- factor(panel$X1, levels = Var$Var)
  panel$X2 <- factor(panel$X2, levels = Var.part$Var)
  for (i in 1:nrow(Var.part)) {
    i <- as.character(Var.part[i, "Var"])
    for (j in unlist(strsplit(i, ", "))) panel[which(panel$X1 ==
                                                       j & panel$X2 == i), "X3"] <- "1"
  }
  panel$X4 <- apply(panel, 1, function(x) ifelse(which(levels(panel$X1) ==
                                                         x[1])%%2 == 0, "0", "1"))


  Var.exp1<-left_join(panel,Hier.part,by=c("X1"="Var"))
  Var.exp1$yaxis<-Var.exp1$X1
  Var.exp1<- Var.exp1[order(Var.exp1$Individual,decreasing = TRUE),]
  Var.exp1$X1 <- factor(Var.exp1$X1, levels = unique(Var.exp1$X1))

  Var.exp1$X4 <- apply(Var.exp1, 1, function(x) ifelse(which(levels(Var.exp1$X1) ==
                                                               x[1])%%2 == 0, "0", "1"))

  #Var.exp1$X1<-factor(Var.exp1$X1, levels = unique(Var.exp1$X1)[length(unique(Var.exp1$X1)):1] )
  Hier.part1<-Hier.part
  Hier.part1<-Hier.part1[order(Hier.part1$Individual,decreasing = TRUE),]
  Hier.part1$Var<-factor(Hier.part1$Individual,levels = unique(Hier.part1$Individual))

  p3 <- ggplot(data = Hier.part1,
               aes(x = Var, y = Individual)) +
    geom_text(aes(y = 0,x = Var,
                  label = scales::percent(Individual,
                                          scale = 1,
                                          suffix="",
                                          accuracy = 0.01)),
              family = "serif",
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          text = element_text(family = "serif"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())

  Hier.part2<-tibble::rownames_to_column(Hier.part1,var = "X1")
  Hier.part2$X1<-factor(Hier.part2$X1,levels = unique(Hier.part2$X1)[length(Hier.part2$X1):1])
  p4 <- ggplot(data = Hier.part2,
               aes(x = Var, y = Individual)) +
    geom_text(aes(y = 0,x = Var,
                  label = X1),
              family = "serif",
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          text = element_text(family = "serif"),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())

  Var.exp1$X2 <- factor(Var.exp1$X2, levels = unique(Var.exp1$X2[length(Var.exp1$X2):1]))


  if (length(point.shape)==1) {
    point.shape<-c(point.shape,point.shape)
  } else if (length(point.shape)==2) {
    point.shape<-point.shape
  } else if (length(point.shape)>2) {
    stop("Please provide no more than two shapes !!!")
  }

  point.shape.use<-point.shape
  names(point.shape.use)<-unique(Var.exp1$X3)[order(unique(Var.exp1$X3))]
  p.panel <- ggplot2::ggplot(data = Var.exp1,
                             aes(x = X2, y = X1, color = X3, fill = X4,shape=X3)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::geom_point(size = pch.size) +
    ggplot2::scale_shape_manual(values = point.shape.use)+
    ggplot2::scale_fill_manual(values = panel.color1,
                               limits = c("0", "1")) +
    ggplot2::scale_color_manual(values = panel.color2,
                                limits = c("0", "1")) +
    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(),
                   axis.text = element_text(color = "black",
                                            size = title.cex),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::labs(y = NULL, x = NULL)+
    ggplot2::theme(axis.text.y.left = element_blank())



  for (i in levels(Var.exp1$X2)) {
    panel.i <- Var.exp1[which(Var.exp1$X2 == i & Var.exp1$X3 == 1),
    ]
    panel.i <- panel.i[order(panel.i$X1), ]
    if (nrow(panel.i) > 1)
      p.panel <- p.panel + ggplot2::annotate("segment",
                                             x = i, xend = i, y = panel.i[1, 1], yend = panel.i[nrow(panel.i),
                                                                                                1], color = "black", size = line.lwd)
  }


  p.vp<-p.vp+theme(text=element_text(family = "serif"))
  p.panel<-p.panel+theme(text=element_text(family = "serif"))


  if (plot.hp) {
    p.vp<-p.vp +
      annotate("text",  x = 1, y = min(Var.part$Fractions)*1.05, family="serif",
               label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))

    xp<-p.panel +p3 +patchwork::plot_layout(widths=c(9.5,0.5),byrow=TRUE)
    p.vp +xp+  + patchwork::plot_layout(heights = c(3.5,6.5),byrow=TRUE)
  } else {
    p.vp<-p.vp+
      annotate("text",  x = 11,
               y = max(Var.part$Fractions)*1.05,
               family="serif",
               label = paste("Constrained:", Constrained, "%", "  ", "Residual", 100 - Constrained, "%"))

    addp<-ggplot(data = Hier.part2,
                 aes(x = Var, y = Individual)) +
      geom_text(aes(y = 0,x = Var,
                    label = X1),color="white",
                family = "serif",
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            text = element_blank(),axis.ticks = element_blank())

    (p3+p.panel+patchwork::plot_layout(widths = c(0.5,9.5)))/
    (addp+p.vp+patchwork::plot_layout(widths = c(0.5,9.5)))+
    patchwork::plot_layout(heights = c(3.5,6.5))

  }


}


"plotMicrochatRDA" <- function(microchatRDAobj,
                               point.shape=point.shape,
                               axis.layout=c("left-bottom","left-top",
                                             "right-bottom","right-top"),
                               nVar=30,
                               order.part=c('degree','effect'),
                               main.barcolor=colorCustom(2,pal = "gygn")[1:2],
                               panel.color1=c("white", "gray90"),
                               panel.color2=c("gray80", "black"),
                               export_path="microbial canonical analysis") {
  axis.layout<-match.arg(axis.layout)
  order.part<-match.arg(order.part)

  dir.create(export_path, recursive = TRUE)
  if (class(microchatRDAobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  mod1<-microchatRDAobj
  if (axis.layout=="left-bottom") {
    tp<-rda1Visual(mod1,
                   point.shape=point.shape,
                   main.barcolor=main.barcolor,
                   panel.color1=panel.color1,
                   panel.color2=panel.color2,
                   plot.hp = TRUE,
                   order.part = order.part,
                   nVar = nVar)+
      theme(text = element_text(family = "serif"))
  }

  if (axis.layout=="right-bottom") {
    tp<-rda3Visual(mod1,
                   point.shape=point.shape,
                   main.barcolor=main.barcolor,
                   panel.color1=panel.color1,
                   panel.color2=panel.color2,
                   plot.hp = TRUE,
                   order.part = order.part,
                   nVar = nVar)+
      theme(text = element_text(family = "serif"))
  }

  if (axis.layout=="left-top") {
    tp<-rda2Visual(mod1,
                   point.shape=point.shape,
                   main.barcolor=main.barcolor,
                   panel.color1=panel.color1,
                   panel.color2=panel.color2,
                   plot.hp = FALSE,
                   order.part = order.part,
                   nVar = nVar)+
      theme(text = element_text(family = "serif"))
  }

  if (axis.layout=="right-top") {
    tp<-rda4Visual(mod1,
                   point.shape=point.shape,
                   main.barcolor=main.barcolor,
                   panel.color1=panel.color1,
                   panel.color2=panel.color2,
                   plot.hp = FALSE,
                   order.part = order.part,
                   nVar = nVar)+
      theme(text = element_text(family = "serif"))
  }

  ggsave(paste(export_path,"/canonical analysis (",axis.layout,").pdf",sep = ""),tp)
  return(tp)
}

"plotMicrochatComplexRDA" <- function(microchatRDAobj,microchatRDAobj1,
                                      point.shape=point.shape,
                                      layout=c("set1","set2","set3","set4"),
                                      nVar=30,
                                      main.barcolor=colorCustom(2,pal = "gygn")[1:2],
                                      panel.color1=c("white", "gray90"),
                                      panel.color2=c("gray80", "black"),
                                      export_path="microbial canonical analysis") {

  layout<-match.arg(layout)
  dir.create(export_path, recursive = TRUE)
  if (class(microchatRDAobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  if (class(microchatRDAobj1)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  mod1<-microchatRDAobj
  mod2<-microchatRDAobj1

  tp1<-rda1Visual(mod1,
                  point.shape=point.shape,
                  main.barcolor=main.barcolor,
                  panel.color1=panel.color1,
                  panel.color2=panel.color2,
                  plot.hp = TRUE,
                  order.part = "degree",
                  nVar = nVar)+
    theme(text = element_text(family = "serif"))



  tp3<-rda3Visual(mod1,point.shape=point.shape,
                  main.barcolor=main.barcolor,
                  panel.color1=panel.color1,
                  panel.color2=panel.color2,
                  plot.hp = TRUE,
                  order.part = "effect",
                  nVar = nVar)+
    theme(text = element_text(family = "serif"))



  tp2<-rda2Visual(mod2,point.shape=point.shape,
                  main.barcolor=main.barcolor,
                  panel.color1=panel.color1,
                  panel.color2=panel.color2,
                  plot.hp = FALSE,
                  order.part = "degree",
                  nVar = nVar)+
    theme(text = element_text(family = "serif"))



  tp4<-rda4Visual(mod2,point.shape=point.shape,
                  main.barcolor=main.barcolor,
                  panel.color1=panel.color1,
                  panel.color2=panel.color2,
                  plot.hp = FALSE,
                  order.part = "effect",
                  nVar = nVar)+
    theme(text = element_text(family = "serif"))

  if (layout=="set1")
    tp<-cowplot::plot_grid(tp3,tp1,
                           tp4,tp2,ncol = 2)

  if (layout=="set2")
    tp<-cowplot::plot_grid(tp1,tp3,
                           tp2,tp4,ncol = 2)

  if (layout=="set3")
    tp<-cowplot::plot_grid(tp2,tp4,
                           tp1,tp3,ncol = 2)

  if (layout=="set4")
    tp<-cowplot::plot_grid(tp4,tp2,
                           tp3,tp1,ncol = 2)

  ggsave(paste(export_path,"/canonical analysis ( ComplexRDA",").pdf",
               sep = ""),tp)
  return(tp)
}

"calcMicrochatRDA" <- function(submchat,
                               distance=c("bray","jaccard"),
                               method='RDA',
                               sta.method= 'hellinger',
                               paramfile.select=c("gut_gene","gut_enzyme")) {
  distance<-match.arg(distance)
  if (class(submchat) !=  "microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  if (is.null(submchat$param_table)) {stop("\n","Please provide the parametric table !!!")}

  if (distance=="bray") {
    if (class(submchat$param_table)!="list") param<-submchat$param_table
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
      param<-param.select
      }
    if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
    if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

    abun<-submchat$otu_table
    data_abun<-abun%>%t()%>%data.frame()
    data_env<-param

    mite.hel <- vegan::decostand(data_abun, method = sta.method)
    mod1 <- rdacca.hp::rdacca.hp(mite.hel, data_env, method = method, var.part = TRUE, type = 'adjR2', scale = FALSE)
  }

  if (distance=="jaccard") {
    if (class(submchat$param_table)!="list") param<-submchat$param_table
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
      param<-param.select
    }

    microchatcomobj<-calcMicrochatcom(submchat,
                                      taxa="Class",
                                      taxa_num=10,
                                      export_path="microbial composition")
    abun1<-microchatcomobj$plot_data
    abun1<-subset(abun1,select = c(sample,Taxo,value))
    abun2<-tidyr::spread(abun1,key=Taxo,value=value)%>%
      tibble::column_to_rownames(var = "sample")

    data_abun<-abun2%>%data.frame()
    data_env<-param

    mite.hel <- vegan::decostand(data_abun, method = sta.method)
    mod1 <- rdacca.hp::rdacca.hp(mite.hel, data_env, method = method, var.part = TRUE, type = 'adjR2', scale = FALSE)
  }
  microchatRDAobj <-mod1
  class(microchatRDAobj) <- "microchat"

  return(invisible(microchatRDAobj))
}


