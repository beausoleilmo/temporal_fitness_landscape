loading <- function(rdata_file)
{
  load(rdata_file, envir = globalenv(),verbose=TRUE)
}

load.data <- function(model ="_pc2", 
                      path_model = NULL,
                      path_fn = "~/Dropbox/finch_recap/v2/saved/",
                      spline.logic = FALSE, 
                      yr.l = 2003:2014, 
                      sp.l =  c('fortis')#, 'fuliginosa', 'scandens', 'magnirostris')
                      ) {
sp.list <- sp.l
yr.list <- yr.l
splines <- spline.logic
setwd('~/Dropbox/finch_recap/v2')
source('src/initialize.R')
## ----------------------------------------------------------------------
loading('data/data.RData')
## ----------------------------------------------------------------------
## prepare data
mod <<- model

spl <<- ifelse(splines, '_splines', '')
if(is.null(path_model)){
  if(length(sp.list)==1) {
    source(sprintf('src/model_ss%s.R', mod))
    fn <<- sprintf('out_ss%s%s.RData', mod, spl)
  }
  if(length(sp.list)>1) {
    source(sprintf('src/model_ms%s.R', mod))
    fn <<- sprintf('out_ms%s%s.RData', mod, spl)
  }
  loading(file.path('saved', fn))
} else {
  # setwd(path_model)
  if(length(sp.list)==1) {
    source(file.path(path_model,sprintf('model_ss%s.R', mod)))
    fn <<- sprintf('out_ss%s%s.RData', mod, spl)
  }
  if(length(sp.list)>1) {
    source(file.path(path_model,sprintf('model_ms%s.R', mod)))
    fn <<- sprintf('out_ms%s%s.RData', mod, spl)
  } 
  path_fn = path_fn
  loading(file.path(path_fn, fn))
}

cols <<- c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')

summ <<- model.summary$BUGSoutput$summary

return(list(fn = fn, 
            # summ = summ,
            cols = cols))
}

get.rows <- function(ss, nn)
  sapply(1:nn, function(x) sprintf(ss, x))

get.rowsz <- function(ss, 
                      nind = my.data$nind, 
                      nyear = my.data$nyear) {
sapply(1:nyear, function(yr)
  sapply(1:nind, 
         function(ind) sprintf(ss, 
                               yr, ind)))
}

plot.fitness = function(xlab.p = 'pc2 value',
                        ylim.phi=c(-2,5),
                        pdf.print = TRUE) {
  if(pdf.print){
  pdf(paste("~/Dropbox/finch_recap/v2/figures/fit.land/fitness_",
            gsub(".RData","",gsub("out_","",fn)),
            ".pdf",sep = ""))}
  plot(NA, 
       main = paste("fitness_",gsub(".RData","",gsub("out_","",fn)), sep = ""),
       xlim = range(my.data$pc2), 
       ylim = ylim.phi,
       xlab = xlab.p, 
       ylab = 'Survival probability', 
       las=1)
  
  add.yr <- function(yr) {
    rows <- match(sapply(1:my.data$nind,
                         function(x) sprintf('e.phi[%d,%d]', x, yr)),
                  rownames(summ))
    points(x = my.data$pc2[my.data$X[,yr]==1],
           y = summ[rows[my.data$X[,yr]==1],'mean'],
           pch = 16, 
           # ylim = c(0,1), 
           cex=0.5)
    vals <- sapply(1:my.data$nind, 
                   calc.exp, yr=yr)
    lines(x = my.data$pc2[order(my.data$pc2)], 
          y = vals[order(my.data$pc2)],
          col = 'red')
  }
  invisible(sapply(1:my.data$nyear, add.yr))
  if(pdf.print){dev.off()}
  # return(vals = vals)
}






plot.fitness.3d = function(xlab.p = 'pc1 value',
                           ylab.p = 'pc2 value',
                           yr = 1,
                        # ylim.phi=c(-2,5),
                        pdf.print = TRUE) {
  library(rgl)
  if(pdf.print){
    pdf(paste("~/Dropbox/finch_recap/v2/figures/fit.land/fitness_",
              gsub(".RData","",gsub("out_","",fn)),
              ".pdf",sep = ""))}
  plot3d(x = my.data$pc1[my.data$X[,yr]==1],
         y = my.data$pc2[my.data$X[,yr]==1],
         z =  expit(summ[sprintf('e.phi[%d,%d]', 1:my.data$nind,yr),"mean"][my.data$X[,yr]==1]),
       main = paste("fitness_",gsub(".RData","",gsub("out_","",fn)), sep = ""),
       # xlim = range(my.data$pc1), 
       # ylim = range(my.data$pc2),
       xlab = xlab.p, 
       ylab = ylab.p#, 
       # las=1
       )
  if(pdf.print){dev.off()}
  # return(vals = vals)
  add.surface(yr=yr)
}

add.surface <- function(sp.l =  c('fortis'),#, 'fuliginosa', 'scandens', 'magnirostris')
                        yr = 1) {
  trait.phi = cbind(my.data$pc1[my.data$X[,yr]==1],
                    my.data$pc2[my.data$X[,yr]==1],
                    expit(summ[sprintf('e.phi[%d,%d]', 1:my.data$nind,yr),"mean"][my.data$X[,yr]==1]))
  colnames(trait.phi) = c("pc1","pc2","phi")
  trait.phi = as.data.frame(trait.phi)
  
  for(j in 1:length(sp.l)){
    # df3d.sp = df3d[df3d$sp == sp[j],]
    if(nrow(trait.phi) == 1){next} else{
      
      library(akima)
      s = interp(x = trait.phi$pc1,
                 y = trait.phi$pc2,
                 z = trait.phi$phi,
                 duplicate="strip",
                 linear = TRUE,
                 nx = 25,
                 ny = 25)
      #Create a function to generate a continuous color palette
      rbPal <- colorRampPalette(c('red','yellow'))
      #This adds a column of color values
      # based on the y values
      nb.div = 100
      data.col = as.data.frame(matrix(as.factor(cut(s$z,ordered_result = T,
                                                    include.lowest = TRUE,
                                                    right = TRUE, 
                                                    breaks = nb.div)),
                                      dim(s$z)[1],
                                      dim(s$z)[2],byrow = FALSE))
      
      range = levels(cut(s$z,ordered_result = T,
                         include.lowest = TRUE,
                         right = TRUE, 
                         breaks = nb.div))
      library(plyr)
      for(i in 1:ncol(data.col)){
        data.col[,i] <- mapvalues(data.col[,i], 
                                  from=range, 
                                  to=rbPal(nb.div), 
                                  warn_missing = FALSE)
        
      }
      
      
      # col.index = matrix(as.numeric(unlist(data.col)),
      #                    dim(s$z)[1],dim(s$z)[2], byrow = FALSE)
      # 
      # Col <- rbPal(nb.div)[col.index]
      # col= matrix(Col, dim(s$z)[1], dim(s$z)[2])
      ###
      
      # bgplot3d( suppressWarnings ( image.plot( legend.only=TRUE, legend.args=list(text='legend'), zlim=Hlim,col=col) )  )#legend
      surface3d(x = s$x,y = s$y,z = s$z,
                color = as.matrix(data.col), 
                alpha = .8)
    }
  }
}
