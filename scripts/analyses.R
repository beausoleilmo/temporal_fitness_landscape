analyse <- function(d,
                    ni,
                    nt,
                    nb,
                    nc,
                    save.path, 
                    init=NULL) {

  cat(sprintf('Started at %s.\n', format(Sys.time(), "%X")))
  start.time <- Sys.time()
  ptm <- proc.time() # timer for the code
  
  my.data <- list(X      = d$X,
                  nind   = nrow(d$X),
                  nyear  = ncol(d$X),
                  nyr  = ncol(d$X),
                  ngr    = length(unique(d$ind.vars$gr.allsp)),
                  gr.sp  = d$ind.vars$gr.allsp,
                  band   = d$ind.vars$band,
                  mbl    = d$ind.vars$mbl,
                  mbw    = d$ind.vars$mbw,
                  mbd    = d$ind.vars$mbd,
                  mbl.non.scaled = d$ind.vars$mbl.non.scaled,
                  mbw.non.scaled = d$ind.vars$mbw.non.scaled,
                  mbd.non.scaled = d$ind.vars$mbd.non.scaled,
                  mbl_sq = d$ind.vars$mbl_sq,
                  mbw_sq = d$ind.vars$mbw_sq,
                  mbd_sq = d$ind.vars$mbd_sq,
                  mbl_cu = d$ind.vars$mbl_cu,
                  mbw_cu = d$ind.vars$mbw_cu,
                  mbd_cu = d$ind.vars$mbd_cu,
                  mass   = d$ind.vars$mass,
                  tars   = d$ind.vars$tars,
                  wingC  = d$ind.vars$wingC,
                  pc1.bod= d$ind.vars$pc1.bod,
                  pc2.bod= d$ind.vars$pc2.bod,
                  pc1    = d$ind.vars$pc1_lin,
                  pc2    = d$ind.vars$pc2_lin,
                  pc1_sq = d$ind.vars$pc1_sq,
                  pc2_sq = d$ind.vars$pc2_sq,
                  pc1_cu = d$ind.vars$pc1_cu,
                  pc2_cu = d$ind.vars$pc2_cu,
                  first  = d$first, ## first is the visit that each most
                                    ## was first captured on
                  sp = as.numeric(d$ind.vars$sp), ## Levels:
                                                  ## 1=fortis
                                                  ## 2=fuliginosa
                                                  ## 3=magnirostris
                                                  ## 4=scandens
                  age = as.numeric(d$ind.vars$age), ## Levels: 
                                                  ## 1 = f 
                                                  ## 2 = j 
                                                  ## 3 = m
                  nage = length(levels(d$ind.vars$age)),
                  nsp     = length(levels(as.factor(d$ind.vars$sp))),
                  effort  = d$yr.vars$eff.cap.hour.std,
                  ## abundance.year = d$yr.vars$abundance.year,
                  ## abund.sp.year  = d$yr.vars$abund.sp.year,
                  precipitation.year = d$yr.vars$precipitation.year,
                  ## bin.rainfall = d$yr.vars$bin.rainfall,
                  sum.rainfall = d$yr.vars$sum.rainfall,
                  splines = as.numeric(splines),
                  Xinit = d$splines.vars$Xinit,
                  Xfix = d$splines.vars$Xfix, # define the fixed
                                        # effect matrices
                                        # X variables with a column of 1s
                  Xsigma = d$splines.vars$Xsigma, 
                  Z = d$splines.vars$Z, # Z matrix from spline analysis 
                  nknotsb = d$splines.vars$nknotsb, # Number of knots 
                  inits.spline = d$splines.vars$inits.spline, #first and second list of inits
                  # If you want to redifine center per year, JAGS doesn't like to have these define as 1 space to put values in. So I have to let these variables go. If you want to use the center of the phenotypic distribution, you can use these values again  
                  # center1 = d$ind.vars$center1[1],
                  # center2 = d$ind.vars$center2[1],
                  indicator =d$ind.vars$indicator) 

  z.init <- d$X
  z.init[z.init==0] <- 1
  
  ## create inits for JAGS
  my.inits <- function() {
    c(list(z=z.init), init)
  }
  
  ## specify the parameters to be traced by JAGS (this function is
  ## defined in src/model.R) # at the end completly 
  my.params <- get.params()
  
  ## package the data
  dd <- list(data=my.data,
             inits=my.inits,
             params=my.params)

  ## run the model (see src/model.R)
  model.out <- run.model(d=dd,
                         ni=ni,
                         nt=nt,
                         nb=nb,
                         nc=nc)
  
  ## save the output, along with the data used to create it
  model.summary <- model.out #$BUGSoutput$summary
  
  cat(sprintf('Completed at %s.\n', format(Sys.time(), "%X")))
  format(Sys.time(), "%X")
  endoftime <- proc.time() - ptm # end of the timer # look at elapsed!
  endoftime.min <- endoftime[3]/60 #endoftime/60 -> voir elapsed = min
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(sprintf('Run took %f minutes\n', endoftime.min))
  
  save(my.data, 
       model.summary,
       file=save.path,
       duration=endoftime.min)

  NULL
}
