# Custom PCA plot ---------------------------------------------------------

custom_pca <- function(pca = pca, 
                       centered = FALSE, # If you want to separate a points axis and a vector axis 
                       main = 'PCA', # Name of the graph
                       scaling = 2,
                       choices = c(1,2),
                       zoom = 1,
                       square = FALSE,
                       square.size = 1,
                       threshold = 0, 
                       pch = 4,
                       # This is in relation to the vectors 
                       col.arrow = "grey",
                       col.label = "red",
                       size.label = .7,
                       
                       col.sites = "black",
                       shape.points = c(17,16,15,19),
                       size.points = 0.5,
                       size.sites = .7,
                       label.the.arrows = TRUE,
                       vector.names = NULL,
                       
                       label.the.points = FALSE,
                       label.pt = NULL,
                       position.text.arrow =1,
                       # arrow.centered = FALSE,
                       vector.to.colour.col = 1, #c(1:dim(summary(pca)$sites)[1]),
                       vector.to.colour.pch = 1, #c(1:dim(summary(pca)$sites)[1]),
                       automatic.col = TRUE,
                       colour.vector = c("#FF3030","#9ACD31","#1D90FF","#FF8001"),
                       ordiellipse = FALSE,
                       ordiellipse.vector = NULL,
                       rev.x = FALSE, # not done yet
                       rev.y = FALSE, # not done yet
                       transparent.back = FALSE,
                       arrows.different.axis = FALSE,
                       length.arrow = 1, 
                       pca.circle = FALSE,
                       label.figure = NULL,
                       small = NULL,
                       png = NULL,
                       png.w = 9,
                       png.h = 9,
                       res.dpi = 300,
                       pdf = NULL) { 
  if (!is.null(png)) {
    png(file = png, 
        width = png.w, height = png.h, res = res.dpi, units = "in")
  }
  if (!is.null(pdf)) {
    pdf(file = pdf, 
        width = 8.5, height = 11)
  }
  
  summ.pca = summary(pca, scaling = scaling) # Valeurs propres, position des espèces et des sites
  pca.scores1 = scores(pca, scaling = scaling, choices = choices)
  prop.exp = round(summ.pca$cont$importance["Proportion Explained",] *100, 2)
  pca.axis.name = c(paste("PC1 ", "(", prop.exp[1], "%", ")", sep = ""),
                    paste("PC2 ", "(", prop.exp[2], "%", ")", sep = ""),
                    paste("PC3 ", "(", prop.exp[3], "%", ")", sep = ""))[choices]
  zoom = zoom
  threshold = threshold # this threshold is to put in the graph only the plants that are far away from the center 
  # (easier to see on the graph), so the one that contributs the most to the variation 
  
  "pcacircle" <- function (pca) 
  {
    # Draws a circle of equilibrium contribution on a PCA plot 
    # generated from a vegan analysis.
    # vegan uses special constants for its outputs, hence 
    # the 'const' value below.
    
    eigenv <- pca$CA$eig
    p <- length(eigenv)
    n <- nrow(pca$CA$u)
    tot <- sum(eigenv)
    const <- ((n - 1) * tot)^0.25
    radius <- (2/p)^0.5
    radius <- radius * const
    symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)
  }

  line2user <- function(line, side) {
    lh <- par('cin')[2] * par('cex') * par('lheight') *.5
    x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
           `1` = grconvertY(-line * y_off, 'npc', 'user'),
           `2` = grconvertX(-line * x_off, 'npc', 'user'),
           `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
           `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
           stop("Side must be 1, 2, 3, or 4", call.=FALSE))
  }
  
  addfiglab <- function(lab, xl = par()$mar[2], yl = par()$mar[3]) {
    
    text(x = line2user(xl, 2), y = line2user(yl, 3), 
         lab, xpd = NA, font = 2, cex = 1.5, adj = c(0, 1))
    
  }
  
  # costum alignment of the axes
  new_lim <- function(a, type = 1) {
    newdata_ratio <-  NULL
    i <- type * 2 - 1
    old_lim <- par("usr")[i:(i+1)] + c(diff(par("usr")[i:(i+1)]) * 0.04 / 1.08, 
                                       diff(par("usr")[i:(i+1)]) * -0.04 / 1.08)
    old_ratio <- old_lim[1] / old_lim[2]
    newdata_ratio <- if (max(a) <= 0) -1.0e+6 else min(a) / max(a)
    if (old_ratio >= newdata_ratio ) {
      new_min <- min(a)
      new_max <- min(a) / old_ratio
    } else {
      new_min <- max(a) * old_ratio
      new_max <- max(a)
    }
    c(new_min, new_max)
  }
  
  new_xylim <- function(x2, y2) {
    
    new_lim2 <- function(a, type = 1, ratio = FALSE) {
      new_ratio <-  NULL
      i <- type * 2 - 1
      old_lim <- par("usr")[i:(i+1)] + c(diff(par("usr")[i:(i+1)]) * 0.04 / 1.08, 
                                         diff(par("usr")[i:(i+1)]) * -0.04 / 1.08)
      old_ratio <- diff(old_lim) / (-1 * old_lim[1])
      new_ratio <- if (min(a) >= 0) 1.0e+6 else diff(range(a)) / (-1 * min(a))
      if (old_ratio >= new_ratio) {
        new_min <- min(a)
        new_max <- abs(min(a)) * (old_ratio - 1)
      } else {
        new_min <- abs(max(a)) / (old_ratio - 1) * -1
        new_max <- max(a)
      }
      if (ratio) return(old_ratio) else return(c(new_min, new_max))
    }
    
    temp_x <- new_lim2(x2, type=1)  # calculation of needing range
    temp_y <- new_lim2(y2, type=2)
    
    if (par("pin")[1] / diff(temp_x) >= par("pin")[2] / diff(temp_y)) {
      # logical judgment of which range is expanded under the same scale
      new_ylim <- temp_y
      new_xdiff <- par("pin")[1] * diff(range(new_ylim)) / par("pin")[2] # calculation of diff(xrange) under the same scale
      old_xratio <- new_lim2(x2, type = 1, ratio = TRUE)
      new_xlim <- c(-1 * new_xdiff / old_xratio, new_xdiff - new_xdiff / old_xratio) # decision of xlim by first graph's zero position and new_xdiff
    } else {
      new_xlim <- temp_x
      new_ydiff <- par("pin")[2] * diff(range(new_xlim)) / par("pin")[1]
      old_yratio <- new_lim2(y2, type=2, ratio = TRUE)
      new_ylim <- c(-1 * new_ydiff / old_yratio, new_ydiff - new_ydiff / old_yratio)
    }
    return(list(new_xlim, new_ylim))
  }
  

# Draw the PCA with plot.cca
  
myccaplot <- function (x, choices = c(1, 2), display = c("sp", "wa", "cn"), 
                         scaling = "species", type, xlim, ylim, const, correlation = FALSE, 
                         hill = FALSE, ...) {
    TYPES <- c("text", "points", "none")
    g <- scores(x, choices, display, scaling, const, correlation = correlation, 
                hill = hill)
    if (length(g) == 0 || all(is.na(g))) 
      stop("nothing to plot: requested scores do not exist")
    if (!is.list(g)) 
      g <- list(default = g)
    for (i in seq_along(g)) {
      if (length(dim(g[[i]])) > 1) 
        rownames(g[[i]]) <- rownames(g[[i]], do.NULL = FALSE, 
                                     prefix = substr(names(g)[i], 1, 3))
    }
    if (!is.null(g$centroids)) {
      if (is.null(g$biplot)) 
        g$biplot <- scores(x, choices, "bp", scaling)
      if (!is.na(g$centroids)[1]) {
        bipnam <- rownames(g$biplot)
        cntnam <- rownames(g$centroids)
        g$biplot <- g$biplot[!(bipnam %in% cntnam), , drop = FALSE]
        if (nrow(g$biplot) == 0) 
          g$biplot <- NULL
      }
    }
    if (missing(type)) {
      nitlimit <- 80
      nit <- max(nrow(g$spe), nrow(g$sit), nrow(g$con), nrow(g$def))
      if (nit > nitlimit) 
        type <- "points"
      else type <- "text"
    }
    else type <- match.arg(type, TYPES)
    if (length(choices) == 1) {
      if (length(g) == 1) 
        pl <- linestack(g[[1]], ...)
      else {
        hasSpec <- names(g)[1] == "species"
        ylim <- range(c(g[[1]], g[[2]]), na.rm = TRUE)
        pl <- linestack(g[[1]], ylim = ylim, side = ifelse(hasSpec, 
                                                           "left", "right"), ...)
        linestack(g[[2]], ylim = ylim, side = ifelse(hasSpec, 
                                                     "right", "left"), add = TRUE, ...)
      }
      return(invisible(pl))
    }
    if (missing(xlim)) {
      xlim <- range(g$species[, 1], g$sites[, 1], g$constraints[, 
                                                                1], g$biplot[, 1], if (length(g$centroids) > 0 && 
                                                                                       is.na(g$centroids)) NA else g$centroids[, 1], g$default[, 
                                                                                                                                               1], na.rm = TRUE)
    }
    if (!any(is.finite(xlim))) 
      stop("no finite scores to plot")
    if (missing(ylim)) {
      ylim <- range(g$species[, 2], g$sites[, 2], g$constraints[, 
                                                                2], g$biplot[, 2], if (length(g$centroids) > 0 && 
                                                                                       is.na(g$centroids)) NA else g$centroids[, 2], g$default[, 
                                                                                                                                               2], na.rm = TRUE)
    }
    if(rev.x) {
      xlim=rev(xlim)
      }
    if(rev.y) {
      ylim=rev(ylim)
    }
    plot(g[[1]], xlim = xlim, ylim = ylim, type = "n", 
         asp = 1, # I had to remove this so that new_lim work.
         ...)
    # abline(h = 0, lty = 3) # removed
    # abline(v = 0, lty = 3) # removed
    if (!is.null(g$species)) {
      if (type == "text") 
        text(g$species, rownames(g$species), col = "red", 
             cex = 0.7)
      else if (type == "points") 
        points(g$species, pch = "+", col = "red", cex = 0.7)
    }
    if (!is.null(g$sites)) {
      if (type == "text") 
        text(g$sites, rownames(g$sites), cex = 0.7)
      else if (type == "points") 
        points(g$sites, pch = 1, cex = 0.7)
    }
    if (!is.null(g$constraints)) {
      if (type == "text") 
        text(g$constraints, rownames(g$constraints), cex = 0.7, 
             col = "darkgreen")
      else if (type == "points") 
        points(g$constraints, pch = 2, cex = 0.7, col = "darkgreen")
    }
    if (!is.null(g$biplot) && nrow(g$biplot) > 0 && type != "none") {
      if (length(display) > 1) {
        mul <- ordiArrowMul(g$biplot)
      }
      else mul <- 1
      attr(g$biplot, "arrow.mul") <- mul
      arrows(0, 0, mul * g$biplot[, 1], mul * g$biplot[, 2], 
             length = 0.05, col = "blue")
      biplabs <- ordiArrowTextXY(mul * g$biplot, rownames(g$biplot))
      text(biplabs, rownames(g$biplot), col = "blue")
      axis(3, at = c(-mul, 0, mul), labels = rep("", 3), col = "blue")
      axis(4, at = c(-mul, 0, mul), labels = c(-1, 0, 1), col = "blue")
    }
    if (!is.null(g$centroids) && !is.na(g$centroids) && type != 
        "none") {
      if (type == "text") 
        text(g$centroids, rownames(g$centroids), col = "blue")
      else if (type == "points") 
        points(g$centroids, pch = "x", col = "blue")
    }
    if (!is.null(g$default) && type != "none") {
      if (type == "text") 
        text(g$default, rownames(g$default), cex = 0.7)
      else if (type == "points") 
        points(g$default, pch = 1, cex = 0.7)
    }
    class(g) <- "ordiplot"
    invisible(g)
  }
 
# I don't remember why I have this 
  align <- function(x) {
      c(-max(c(-min(x),max(x))),max(c(-min(x),max(x))))
    }

# My costum function 
  if (transparent.back) {
    par(bg=NA)  
  }
  
  if (centered) {
    par(mar = c(6,6,6,6))
    x0 = c(min(pca.scores1$sites[,1])*zoom,max(pca.scores1$sites[,1])*zoom)
    # pmx = c(-max(abs(x0)),max(abs(x0)))
    y0 = c(min(pca.scores1$sites[,2])*zoom,max(pca.scores1$sites[,2])*zoom)
    # pmy = c(-max(abs(y0)),max(abs(y0)))
    myccaplot(pca,
         # choices = choices,
         # type = c("text",
         # "points"),
         type = c("n"),
         main = main,
         scaling = scaling,
         choices = choices,
         bty = "o",
         xlab = pca.axis.name[1], ylab = pca.axis.name[2],
         xlim = x0, # pmx,
         ylim = y0, # pmy,
         # xlim = c(min(pca.scores1$species[,1],pca.scores1$sites[,1])*zoom,
         #          max(pca.scores1$species[,1],pca.scores1$sites[,1])*zoom),
         # ylim = c(min(pca.scores1$species[,2],pca.scores1$sites[,2])*zoom,
         #          max(pca.scores1$species[,2],pca.scores1$sites[,2])*zoom),
         pch=4,
         axes = FALSE)
    axis(1)  ### (2)
    axis(2)  ### (2)
    box()    ### (2)
    
    points(pca.scores1$sites,
           col=colour.vector[(as.numeric(vector.to.colour.col))], # c("grey","blue") # palette()
           pch=shape.points[(as.numeric(vector.to.colour.pch))], # red = AB, green = EG, light blue = magnirostris, purple = scandens
	   bg= colour.vector[(as.numeric(vector.to.colour.col))],
           cex = size.points)
    if (label.the.points) {
      text(x = summ.pca$sites[,1],
           y = summ.pca$sites[,2],
           labels = label.pt, #row.names(summ.pca$sites[,]), 
           cex = size.sites, 
           col = col.sites, 
           pos = 1)    
    }
    if (ordiellipse) {
      ordiellipse(pca.scores1,ordiellipse.vector,conf=0.95)
    }
      par(new=TRUE)
      x1 = c(min(pca.scores1$species[,1]),max(pca.scores1$species[,1]))
      y1 = c(min(pca.scores1$species[,2]),max(pca.scores1$species[,2]))
      plot(x = x1,
           y = y1,
           xlim = new_xylim(x1, y1)[[1]],
           ylim = new_xylim(x1, y1)[[2]],
           asp = 1, #Not working with new_lim
           xlab = pca.axis.name[1],
           ylab = pca.axis.name[2],
           type = "n",
           main= "",
           axes = F, ann = F,
           bty = "n")
      axis(3, col="red",col.axis="red",
           # las=2, 
           cex.axis = 1)
      axis(4, col="red",col.axis="red",
           # las=2, 
           cex.axis = 1)
      if (pca.circle == TRUE & scaling == 1){
        pcacircle(pca)
      }
      
      # mtext("Arrow PC1", side = 3, 
      #       line = 2.5, # if you increase line, you make the title farther
      #       col = "red")
      # mtext("Arrow PC2", side = 4, 
      #       line = 2.5, 
      #       col = "red")
      abline(h=0,v=0, lty = 3)
      
      if (rev.x) {
        pca.scores1$species[,1] = -1*(pca.scores1$species[,1]*length.arrow)
      }
      if (rev.y) {
        pca.scores1$species[,2] = -1*(pca.scores1$species[,2]*length.arrow)
      }
      arrows(0,0, 
             pca.scores1$species[,1]*length.arrow, 
             pca.scores1$species[,2]*length.arrow, 
             length = .1, 
             code = 2, 
             col = col.arrow)
    
  } else {
    plot(pca,
         choices = choices,
         # type = c("text",
         # "points"),
         type = c("n"),
         xlab = pca.axis.name[1],
         ylab = pca.axis.name[2],
         main= main,
         scaling = scaling,
         # xlim = c(min(pca.scores1$species[,1],pca.scores1$sites[,1])*zoom,
         #          max(pca.scores1$species[,1],pca.scores1$sites[,1])*zoom),
         # ylim = c(min(pca.scores1$species[,2],pca.scores1$sites[,2])*zoom,
         #          max(pca.scores1$species[,2],pca.scores1$sites[,2])*zoom),
         pch=4)
    points(pca.scores1$sites,
           col=colour.vector[(as.numeric(vector.to.colour.col))], # c("grey","blue") # palette()
           pch=shape.points[(as.numeric(vector.to.colour.pch))], # red = AB, green = EG, light blue = magnirostris, purple = scandens
           cex = size.points)
    if (rev.x) {
      summ.pca$species[,1] = rev(summ.pca$species[,1])
    }
    if (rev.y) {
      summ.pca$species[,1] = rev(summ.pca$species[,1])
    }
    arrows(0,0, 
           summ.pca$species[,choices[1]], 
           summ.pca$species[,choices[2]], 
           length = .1, 
           code = 2, 
           col = col.arrow)
    if (pca.circle == TRUE & scaling == 1){
      pcacircle(pca)
    }
    if (ordiellipse) {
      ordiellipse(pca.scores1,vector.to.colour.col,conf=0.95)
    }
  }
  # if (automatic.col) {
  #   points(pca.scores1$sites*inverse, 
  #          col=palette()[(as.numeric(vector.to.colour))], # c("grey","blue") # 
  #          pch=shape.points[(as.numeric(vector.to.colour))],# red = AB, green = EG, light blue = magnirostris, purple = scandens
  #          cex = size.points) 
  #   
  # } else {
    # points(pca.scores1$sites*inverse,
    #        col=colour.vector[(as.numeric(vector.to.colour))], # c("grey","blue") # palette()
    #        pch=shape.points[(as.numeric(vector.to.colour))], # red = AB, green = EG, light blue = magnirostris, purple = scandens
    #        cex = size.points)
  # }
  
  if (rev.x) {
    summ.pca$species[,1] = -1*(summ.pca$species[,1])
  }
  if (rev.y) {
    summ.pca$species[,2] = -1*(summ.pca$species[,2])
  }
  if (!(length(small))) {
    if (label.the.arrows) {
      text(x = summ.pca$species[,choices[1]],
           y = summ.pca$species[,choices[2]],
           labels = row.names(summ.pca$species[,]), 
           cex = size.label, 
           col = col.label)
    } else {
      text(x = summ.pca$species[,choices[1]],
           y = summ.pca$species[,choices[2]],
           labels = vector.names, 
           pos = position.text.arrow,
           cex = size.label, 
           col = col.label)
    }
  } else {
    if (label.the.arrows) {
      text(x = summ.pca$species[-small,choices[1]],
           y = summ.pca$species[-small,choices[2]],
           labels = row.names(summ.pca$species[-small,]), 
           cex = size.label, 
           col = col.label)
    } else {
      text(x = summ.pca$species[,choices[1]],
           y = summ.pca$species[,choices[2]],
           labels = vector.names, 
           pos = position.text.arrow,
           cex = size.label, 
           col = col.label)
    }
  }

  # text(pca, display = "species", choices = c(1, 2), scaling = 2, head.arrow = 0.05, cex = .7, col = "red")
  if (rev.x) {
    summ.pca$sites[,1] = -1*(summ.pca$sites[,1])
  }
  if (rev.y) {
    summ.pca$sites[,2] = -1*(summ.pca$sites[,2])
  }
  if (label.the.points & !centered) {
    text(x = summ.pca$sites[,1],
         y = summ.pca$sites[,2],
         labels = label.pt, #row.names(summ.pca$sites[,]), 
         cex = size.sites, 
         col = col.sites, 
         pos = 1)  
  }
  if (square) {
    abline(h = c(-square.size,square.size), 
           v = c(-square.size,square.size), 
           lty = 3)
  }
  
  
  # if (ordiellipse) {
  #   ordiellipse(pca.scores1,vector.to.colour,conf=0.95)
  # }
  if (!is.null(png)|!is.null(pdf)) {
    dev.off()
    
  }
  if (!(!length(label.figure))) {
  addfiglab(label.figure)
  }
}

# par("pin")[1]/diff(par("usr")[1:2])
# par("pin")[2]/diff(par("usr")[3:4]) # the same
# In this case on my environment (640 x 500), I can resize the 
# width larger and/or the height smaller (change of axis(1) 
# and axis(3)) but can't the width smaller and/or the height larger (change of axis(1) but not axis(3)).




# Example code ------------------------------------------------------------
# Nb of species should match the number of variable in the mu.trait vecotr
# This part is based on the species in 'sp.id'
# mu.trait  =  structure(list(fortis = c(9.826479, 11.638375, 11.091573), fuliginosa = c(6.838141, 
#                                                                                        8.506954, 7.09429), magnirostris = c(13.167606, 14.480028, 15.274225
#                                                                                        ), scandens = c(8.260292, 14.207143, 8.587088)), .Names = c("fortis", 
#                                                                                                                                                    "fuliginosa", "magnirostris", "scandens"), row.names = c("MedianBeakWidth_mean", 
#                                                                                                                                                                                                             "MedianBeakLength_mean", "MedianBeakDepth_mean"), class = "data.frame")
# sigma.trait= structure(list(fortis = c(0.9860465, 0.9385299, 1.7517643), fuliginosa = c(0.1005923, 
#                                                                                         0.2206302, 0.1554465), magnirostris = c(0.9799356, 0.7565252, 
#                                                                                                                                 1.9891933), scandens = c(0.1857845, 0.8934824, 0.2888495)), .Names = c("fortis", 
#                                                                                                                                                                                                        "fuliginosa", "magnirostris", "scandens"), row.names = c("MedianBeakWidth_var", 
#                                                                                                                                                                                                                                                                 "MedianBeakLength_var", "MedianBeakDepth_var"), class = "data.frame")
# 
# 
# nsp=4
# nind = 100
# sp.id = sort(rep(c(1,2,3,4),25))
# ind.trait <- matrix(0, nrow=nind, ncol=dim(mu.trait)[1]) # creating the matrix of trait values 'empty' 
# nb.id.sp.last=1 # counter, maybe not the most elegant way to do this... 
# for (nb.id.sp in 1:nsp) { # For each species
#   nb.id.sp.seq=as.numeric(summary(as.factor(sp.id))[nb.id.sp]) # I want to know how many individuals are there 
#   for (nbtraits in 1:(dim(mu.trait)[1])) { # I'm going to generate new traits baised on the number of traits that I provide in the variable definition
#     traits = rnorm(n=nb.id.sp.seq, # # Generate traits from the matrix of mean and variance of known traits from a normal distribution  
#                    mean = mu.trait[nbtraits,nb.id.sp], # known mean 
#                    sd = sigma.trait[nbtraits,nb.id.sp]) # known  variance
#     ind.trait[nb.id.sp.last:(nb.id.sp.last+nb.id.sp.seq-1),nbtraits] = traits # Adding traits to the matrix
#   } 
#   nb.id.sp.last=nb.id.sp.last+nb.id.sp.seq # This is the counter to add sequencially new traits, again, maybe not the most elegant way... 
# }
# 
# colnames(ind.trait) <- paste(c('trait'),1:dim(mu.trait)[1], sep = '')
# rownames(ind.trait) <- paste('id', 1:nind, sep='.')
# library(vegan)
# pca=rda(ind.trait)
# pca_scores<-scores(pca, choices = c(1,2), scaling = 2) # saving the pca scores
# PC1 = pca_scores$sites[,1] # saving the pca scores 1
# PC2 = pca_scores$sites[,2] # saving the pca scores 2
# ind.trait.PC <- matrix(0, nrow=nind, ncol=2)
# ind.trait.PC[,1] = PC1
# ind.trait.PC[,2] = PC2

# custom_pca(pca = pca,
#            vector.to.colour = sp.id,
#            ordiellipse = TRUE)
