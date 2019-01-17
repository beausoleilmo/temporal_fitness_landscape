  prep.data <- function(sp.keep,
                      yr.keep,
                      site.keep,
                      flip.pc1 = TRUE, # If the PCA axis is not in the same direction as before, you can flip it here 
                      adults.only = FALSE, # If true, it'll keep only females and males (remove juveniles) 
                      center = TRUE, # If you don't want to find the center before hand, put it at FALSE and specify it in the mode
                      gr.sel = NULL, # applicable for fortis only. This will let you select "big" or "small" morph
                      valley.fortis = FALSE,# If this is true, this term will find the valley (between the small morph of fortis and the big morph of fortis). You can run a model on these 2 after that
                      valley.fortis.fuli = FALSE, # If this is true, this term will find the valley (between fuliginosa and the small morph of fortis). You can run a model on these 2 after that
                      peaks.fortis.fuli = FALSE, # If this is true, this term will find the 2 peaks (of fuliginosa and the one for the small morph of fortis). You can run a model on these 2 after that
                      keep.resid = FALSE, # If this is true, it's going to keep the individuals that are putatively resident (==1)
                      splines, # This term is never used here  
                      two.univariate.traits = FALSE, # For spline only. This is computing PC1 and PC2 as 2 independent univariate traits.
                      keep.last.year.data = TRUE, # If you don't want to drop the individuals from the last year 
                      gam.analysis = FALSE,
                      recalculate.pca = FALSE, # Recalculating the PCA for the subsetted data
                      pca.per.sp.only = FALSE, # If TRUE, changes the PCA calculated across all years for all species to PCA calculated across all years for only 
                      K.nots=NULL,
                      traits) {
    # To make some test
    # sp.keep = sp.list <- c('fortis', "fuliginosa")
    # yr.keep = yr.list <- 2003:2014
    # setwd('~/Dropbox/finch_recap/v2')
    # splines <- FALSE
    # source('src/initialize.R')
    # peaks.fortis.fuli = TRUE
    pca.recalc = NULL
    ## ----------------------------------------------------------------------
    ## load data (from .RData)
    # load('data/data.RData', verbose=TRUE)
    # I added other climate data 06 December 2017
    load('data/bird.data.RData', verbose=TRUE)
    
    # I decided to change the species for these guys because they don't cluster properly otherwise 
    bird.data[bird.data$BANDFINAL %in% c("JP1087","SH157_S","UM174","SH123_S"),"Species1"] <- "fortis"
    bird.data[bird.data$BANDFINAL %in% c("JP3846","JP4935","JP737","SH1137"),  "Species1"] <- "fuliginosa"
    bird.data[bird.data$BANDFINAL == "KGSK2147" & bird.data$Species1 == "fuliginosa","Species1"] <- "fortis"
    bird.data[bird.data$BANDFINAL == "JP819","Species1"]<-"fuliginosa"
    bird.data[bird.data$BANDFINAL == "JP3765","Species1"]<-"scandens"
    
    # Finding misidentified fortis/magnirostris    
    # bird.data[bird.data$PC1 >-2 & bird.data$Species1 == "magnirostris","BANDFINAL"]
    # bird.data[bird.data$PC1 >-2 & bird.data$Species1 == "magnirostris","Species1"] <- "fortis"
    # points(-1*bird.data[-1*bird.data$PC1 >-2 & bird.data$Species1 == "magnirostris",c("PC1")],
    #        bird.data[-1*bird.data$PC1 >-2 & bird.data$Species1 == "magnirostris",c("PC2")], 
    #        col = "orange")
    # dput(droplevels(bird.data[-1*bird.data$PC1 >-2 & bird.data$Species1 == "magnirostris",c("BANDFINAL")]))
    list.misidentified.magnirostris.fortis = c("JH1A459", "JP1589", "JP2528", "JP2541", 
                                               "JP3400", "JP3423", "JP3497", "JP3532", 
                                               "JP3598", "JP3607", "JP3626", "JP3656", 
                                               "JP3665", "JP3818", "JP3860", "JP3865", 
                                               "JP3890", "JP4431", "JP4439", "JP4444", 
                                               "JP4453", "JP4517", "JP4548", "JP4744", 
                                               "JP4958", "P739", "PP1A257", "SH1069_L", 
                                               "SH1102_L", "SH1110", "SH1143", "SH1145_L", 
                                               "SH1154_L", "SH192_L", "SH602_L", "SH978_L")
    bird.data[bird.data$BANDFINAL %in% list.misidentified.magnirostris.fortis,"Species1"] <- "fortis"
    
    if(adults.only){
      bird.data = bird.data[bird.data$age %in% c("f","m"),]
    }
    # effort.estimation = apply(bird.data[bird.data$Species1 == "fortis" & 
    #                   bird.data$Site == "Academy Bay",c("y.2003","y.2004","y.2005","y.2006",
    #                                                     "y.2007","y.2008",
    #                                                     "y.2009","y.2010","y.2011","y.2012",
    #                                                     "y.2013","y.2014","y.2015","y.2016", 
    #                                                     "y.2017","y.2018")],2,sum)
    # eff.site = effort[effort$SITE == "AB", "HOUR"]
    # x = effort.estimation[1:length(eff.site)]
    # plot(log(eff.site)~log(x))
    # 
    # new <- data.frame(x = effort.estimation[c(length(eff.site)+1):length(effort.estimation)])
    # new = new
    # pred.lm= predict(lm(log(eff.site)~log(x)), new, se.fit = FALSE)
    # pred.lm = exp(pred.lm)
    # plot(as.vector(c(eff.site,pred.lm))~effort.estimation, xlab = "number birds",ylab = "Net hours")
    # points(effort.estimation[13:16],pred.lm,col ="red")
    # dput(pred.lm)
    eg.net.hours = structure(c(30.2565802161513, 81.4673085182818, 81.4673085182818, 
                     117.180771285717), .Names = c("y.2015", "y.2016", "y.2017", "y.2018"))
    ab.net.hours = structure(c(37.8638868789056, 38.1772072479725, 37.7292504737987, 
                37.5198266878393), .Names = c("y.2015", "y.2016", "y.2017", "y.2018"
                ))
    
  ## reduce to relevant species
  bird.data = droplevels(bird.data[bird.data$Species1 %in% sp.keep,])
  
  bird.data = droplevels(bird.data[bird.data$Site %in% site.keep,])
  
  if (recalculate.pca) {
    library(vegan)
    pca.recalc = vegan::rda(bird.data[,c("MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth")],scale = FALSE)
    pc.scores = scores(pca.recalc,
                       choices = c(1,2), 
                       scaling = 2)
    bird.data$PC1 = pc.scores$sites[,1]
    bird.data$PC2 = pc.scores$sites[,2]
  }
  
  if (pca.per.sp.only) {
    bird.data$PC1 = bird.data$pc1.score.per.sp
    bird.data$PC2 = bird.data$pc2.score.per.sp
    # bird.data$PC1 = bird.data$pc1.score.per.sp.eg
    # bird.data$PC2 = bird.data$pc2.score.per.sp.eg
  }
  
  # If the PCA axis is not in the same direction as before, you can flip it here 
  if(flip.pc1){
    bird.data$PC1 = -1*bird.data$PC1
  }
  
    # Schaub: 
  # It’s possible that there are transients in the data (individuals that are present in the study area only at one occasion) or that survival has an age structure with juvenile survival being lower than adult survival. Classical goodness-of-fit test as can be done in U-CARE are indicative about these potential problems. Models exist that can account for transients or that have an age structure (see e.g. the paper by Roger Pradel (1997) in Biometrics). Deleting the capture histories of all individual captured only once is wrong and would result in a positive bias of the estimate of survival.
  # Pradel 1997 # Transients can be operationally defined as individuals having a zero probability of survival after their initial capture.
  # View(bird.data[bird.data$maxseen>1,])
  # bird.data = bird.data[bird.data$maxseen>1,]
  ch <- colnames(bird.data)[colnames(bird.data) %in%
                            paste('y', yr.keep, sep='.')]
  ## remove birds that were not observed in selected years
  detected <- apply(bird.data[,ch], 1, sum)>0
  bird.data <- bird.data[detected,]

  ## update "first" so that it correctly reflects selected years
  bird.data$first <- apply(bird.data[,ch], 1,
                           function(x) min(which(x==1)))
  
  ## remove birds whose first capture date was in last year (drop a lot of captures in one year)
  if(keep.last.year.data){
    bird.data <- bird.data[bird.data$first<length(ch),]
  }
  if(gam.analysis){
    known.state.cjs <- function(ch) {
      state <- ch
      for(i in 1:dim(ch)[1]){
        n1 <- min(which(ch[i,]==1))
        n2 <- max(which(ch[i,]==1))
        state[i,n1:n2] <- 1
        state[i,n1] <- 1
      }
      state[state==0] <- 0
      state
    }
    bird.data[,ch] = known.state.cjs(bird.data[,ch])
    # This will only work in pairs of years (can't use it if you are doing 2003:2018, only 2003:2004, or any other combination)
    bird.data = bird.data[bird.data[,ch[1]]==1,]
  }
  # Should we keep ONLY resident birds? If keep.resid is true, yes. 
  if(keep.resid){
    bird.data <- bird.data[bird.data$resid.transient ==1,] 
  }
  # length(bird.data$BANDFINAL)
  
  if(length(sp.keep) == 1 & any(sp.keep == "fortis")){ # And that length(sp.keep)==1 (for the moment, it simplifies)
    # Creating 2 groups for fortis 
    quatile.beak.traits = quantile(bird.data$PC1)
    bird.data$gr = if_else(bird.data$PC1 >= quatile.beak.traits[2],"small","big")
    mean.smal  = mean(bird.data[bird.data$gr=="small",  c("PC1")]) # -0.0629293
    mean.big = mean(bird.data[bird.data$gr=="big",c("PC1")]) # -0.3395537
    # plot(bird.data$PC2~bird.data$PC1, col = as.factor(bird.data$gr)); abline(v = pc1.threshold)
    pc1.threshold = as.numeric(quatile.beak.traits[2])
    d1 <- density(bird.data$PC1[bird.data$PC1>pc1.threshold])
    d2 <- density(bird.data$PC1[bird.data$PC1<pc1.threshold])
    center1 <- rep(d1$x[which.max(d1$y)], nrow(bird.data))
    center2 <- rep(d2$x[which.max(d2$y)], nrow(bird.data))
    indicator <- (bird.data$PC1>(center1[1]+center2[1])/2)*1
  } else {
      bird.data$gr = rep(NA, nrow(bird.data))
      center1 = rep(NA, nrow(bird.data))
      center2 = rep(NA, nrow(bird.data))
      indicator = rep(NA, nrow(bird.data))
  }
  if(length(sp.keep) == 4){ # And that length(sp.keep)==1 (for the moment, it simplifies)
    quatile.beak.traits = quantile(bird.data$PC1[bird.data$Species1=="fortis"])
    bird.data$gr.all.sp = as.numeric(bird.data$Species1)+1
    bird.data$gr.all.sp[bird.data$Species1=="fortis"] = if_else(bird.data$PC1[bird.data$Species1=="fortis"] >= quatile.beak.traits[2],
                                                                1,
                                                                2)
  } else {
    bird.data$gr.all.sp = rep(NA, nrow(bird.data))
  }
  
  if(!is.null(gr.sel)){
    bird.data = droplevels(bird.data[bird.data$gr %in% gr.sel,])
    d1 <- density(bird.data$PC1#[bird.data$PC1>pc1.threshold]
                  )
    center1 <- d1$x[which.max(d1$y)]
    center2 = rep(NA, nrow(bird.data))
    indicator <- (bird.data$PC1>(center1))*1
  }
  
  
  if(valley.fortis){
    # The data is extracted like this: 
    bird.data = bird.data[bird.data$PC1 >=mean.big & bird.data$PC1 <=mean.smal,]
    center1 =     center1[1:nrow(bird.data)]
    center2 =     center2[1:nrow(bird.data)]
    indicator = indicator[1:nrow(bird.data)]
    
  }
  
  if(valley.fortis.fuli & all(sp.keep %in% c("fortis", "fuliginosa"))){
    for.data = droplevels(bird.data[bird.data$Species1 %in% "fortis",])
    ful.data = droplevels(bird.data[bird.data$Species1 %in% "fuliginosa",])
    
    quatile.beak.traits=quantile(for.data$PC1)
    for.data$gr = if_else(for.data$PC1 >= quatile.beak.traits[2],"small","big")
    mean.smal  = mean(for.data[for.data$gr=="small",  c("PC1")]) # -0.0629293
    mean.big = mean(for.data[for.data$gr=="big",c("PC1")]) # -0.3395537
    mean.ful = mean(ful.data$PC1)
    # The data is extracted like this: 
    bird.data = bird.data[bird.data$PC1 >=mean.smal & bird.data$PC1 <=mean.ful,] # this will depend on the PCA axis (Be careful of the sign!)
    bird.data$gr = rep(NA, nrow(bird.data))
    center1 = rep(NA, nrow(bird.data))
    center2 = rep(NA, nrow(bird.data))
    indicator = rep(NA, nrow(bird.data))
    
  }
  
  if(peaks.fortis.fuli & all(sp.keep %in% c("fortis", "fuliginosa"))){
    for.data = droplevels(bird.data[bird.data$Species1 %in% "fortis",])
    ful.data = droplevels(bird.data[bird.data$Species1 %in% "fuliginosa",])
    
    quatile.beak.traits=quantile(for.data$PC1)
    for.data$gr = if_else(for.data$PC1 >= quatile.beak.traits[2],"small","big")
    mean.smal  = mean(for.data[for.data$gr=="small",  c("PC1")]) # -0.0629293
    mean.big = mean(for.data[for.data$gr=="big",c("PC1")]) # -0.3395537
    mean.ful = mean(ful.data$PC1)
    # The data is extracted like this: 
    # pc1.threshold = as.numeric(quatile.beak.traits[2])
    d1 <- density(for.data$PC1[for.data$gr == "small"])
    plot(d1)
    d2 <- density(ful.data$PC1[ful.data$Species1=="fuliginosa"])
    plot(d2)
    for.data = for.data[for.data$gr == "small",]
    ful.data$gr = NA
    bird.data = rbind(for.data,ful.data)
    center1 <- d1$x[which.max(d1$y)]
  # -0.07961267
    center2 <- d2$x[which.max(d2$y)]
    indicator <- (bird.data$PC1>(center1+center2)/2)*1
  }
  length(center1)
  length(center2)
  length(indicator)
  length(bird.data$PC2)
  # plot(bird.data$PC2~bird.data$PC1)
  # After Luke verification, I wanted to make sure I wasn't making errors in the script
  # par(mfrow=c(3,1))
  # hist(bird.data$PC1, breaks = 15)
  # mtext(at = c(pc1.threshold+.2),text = round(pc1.threshold,2),col = "red",line = -3)
  # abline(v = pc1.threshold, col='red')
  # abline(v = center1, col='green') # This is correctly identifying the highest mode, though both modes are present on this side of the threshold...
  # hist(bird.data$PC1[indicator==1], add=T, col='grey',breaks = 15)
  # hist(bird.data$PC1[indicator==0], add=T, col='blue',breaks = 15)
  # load('data/bird.data.RData', verbose=TRUE)  
  # bird.data = droplevels(bird.data[bird.data$Species1 %in% sp.keep,])
  # hist(bird.data$PC1, breaks = 15)
  # mtext(at = c(pc1.threshold+.2),text = round(pc1.threshold,2),col = "red",line = -3)
  # abline(v = pc1.threshold, col='red')
  # plot(bird.data$PC2~bird.data$PC1);abline(v = pc1.threshold, col='red')
  # mtext(at = c(pc1.threshold+.2),text = round(pc1.threshold,2),col = "red",line = -3)
  
  
  # center1
  # center2
  ## ## *** Fill Z matrix with known values ***
  ## known.state.cjs <- function(ch) {
  ##   state <- ch
  ##   for(i in 1:dim(ch)[1]){
  ##     n1 <- min(which(ch[i,]==1))
  ##     n2 <- max(which(ch[i,]==1))
  ##     state[i,n1:n2] <- 1 # defining latent state after the first
  ##                         # observation and the last observation
  ##     state[i,n1] <- 1 # replaces all first observation by NA because
  ##                      # the CJS model is *conditional* on first
  ##                      # capture, so the latent state is only defined
  ##                      # after first capture
  ##   }
  ##   state[state==0] <- 0 # each value where an individual is not
  ##                        # observed should take a value of NA because
  ##                        # we don't know if it survived
  ##   state
  ## }
  ## Z <- bird.data[,ch]
  ## Z[,ch] = known.state.cjs(Z)
  
  ## quick scale function
  sc.num <- function(x) as.numeric(scale(x))

  if(!center){
    center1 = NA
    center2 = NA
  }
  ## *** individual-level variables ***
  ind.vars <- data.frame(sp      = bird.data$Species1,
                         gr      = bird.data$gr, 
                         gr.allsp= bird.data$gr.all.sp,
                         mass    = standardize(bird.data$Mass),
                         sex     = bird.data$Sex0,
                         age     = as.factor(bird.data$age),
                         tars    = standardize(bird.data$Tarsus),
                         wingC   = standardize(bird.data$Wing.Chord),
                         mbl     = standardize(bird.data$MedianBeakLength),
                         mbw     = standardize(bird.data$MedianBeakWidth),
                         mbd     = standardize(bird.data$MedianBeakDepth),
                         mbl.non.scaled = bird.data$MedianBeakLength,
                         mbw.non.scaled = bird.data$MedianBeakWidth,
                         mbd.non.scaled = bird.data$MedianBeakDepth,
                         mbl_sq  = poly(standardize(bird.data$MedianBeakLength),3)[,2],
                         mbw_sq  = poly(standardize(bird.data$MedianBeakWidth),3)[,2],
                         mbd_sq  = poly(standardize(bird.data$MedianBeakDepth),3)[,2],
                         mbl_cu  = poly(standardize(bird.data$MedianBeakLength),3)[,3],
                         mbw_cu  = poly(standardize(bird.data$MedianBeakWidth),3)[,3],
                         mbd_cu  = poly(standardize(bird.data$MedianBeakDepth),3)[,3],
                         resid.t = bird.data$resid.transient,
                         pc1.bod = bird.data$PC.body1,
                         pc2.bod = bird.data$PC.body2,
                         pc1     = bird.data$PC1,
                         pc2     = bird.data$PC2,
                         pc1_lin = poly(bird.data$PC1,3)[,1],
                         pc2_lin = poly(bird.data$PC2,3)[,1],
                         pc1_sq  = poly(bird.data$PC1,3)[,2],
                         pc2_sq  = poly(bird.data$PC2,3)[,2],
                         pc1_cu  = poly(bird.data$PC1,3)[,3],
                         pc2_cu  = poly(bird.data$PC2,3)[,3],
                         band    = bird.data$BANDFINAL,
                         center1 = center1,
                         center2 = center2,
                         indicator = indicator)
  std.mean.sd = data.frame(avg.scl.tar = mean(bird.data$Tarsus),
                  avg.scl.mass= mean(bird.data$Mass),
                  avg.scl.wc  = mean(bird.data$Wing.Chord),
                  avg.scl.mbl = mean(bird.data$MedianBeakLength),
                  avg.scl.mbw = mean(bird.data$MedianBeakWidth),
                  avg.scl.mbd = mean(bird.data$MedianBeakDepth),
                  sd.scl.tar  = sd(bird.data$Tarsus),
                  sd.scl.mass = sd(bird.data$Mass),
                  sd.scl.wc   = sd(bird.data$Wing.Chord),
                  sd.scl.mbl  = sd(bird.data$MedianBeakLength),
                  sd.scl.mbw  = sd(bird.data$MedianBeakWidth),
                  sd.scl.mbd  = sd(bird.data$MedianBeakDepth))
                  
# Year-Level variables ----------------------------------------------------
  ## *** year-level variables ***

  ## total sum of precipitation per year
  precip.year <- df.df[,2]
  precip.year = c(precip.year,rep(mean(precip.year),4))
  names(precip.year) <- 2003:2018

  ## precipitation per year between two sampling events (Amount of
  ## rain at a median sampling date)
  precip.btwn.samples <- sum.rainfall[c(-1,-2)]
  precip.btwn.samples = c(precip.btwn.samples,rep(mean(precip.btwn.samples),4))
  names(precip.btwn.samples) <- 2003:2018
  
  # From script Climate_all_time_Santa_cruz.R (On Andrew's demand)
  # > dput(wet.rain.our.study)
  wet.rain.our.study = structure(list(Year.marco = c("2003", "2004", "2005", "2006", 
                                                     "2007", "2008", "2009", "2010", 
                                                     "2011", "2012", "2013", "2014", 
                                                     "2015", 
                                                     "2016", "2017", "2018"), 
                                      sum.prec = c(83.5, 51.1, 123.8, 55.1, 134.8, 664.9, 
                                                   134, 446.2, 611.7, 435.1, 189.2, 176.9, 412.8,
                                                   270.7,270.7,270.7)), # Means
                                 .Names = c("Year.marco", 
                                            "sum.prec"), 
                                 row.names = c(NA, -16L), class = c("tbl_df", "tbl", 
                                                                    "data.frame"))
  wet.rain.our.study = as.data.frame(wet.rain.our.study)
  # > dput(dry.rain.our.study)
  dry.rain.our.study = structure(list(Year.marco = c("2003", "2004", "2005", "2006", 
                                                     "2007", "2008", "2009", "2010", 
                                                     "2011", "2012", "2013", "2014", 
                                                     "2015", 
                                                     "2016", "2017", "2018"), 
                                      sum.prec = c(106.7, 110.1, 62.1, 110, 73.9, 104.1, 85.2, 
                                                   56.6, 67.2, 79.3, 57.3, 71.2, 159.1,
                                                   87.90769,87.90769,87.90769)), 
                                 .Names = c("Year.marco", 
                                            "sum.prec"), 
                                 row.names = c(NA, -16L), class = c("tbl_df", "tbl", 
                                                                    "data.frame"))
  dry.rain.our.study = data.frame(dry.rain.our.study)
  
  
  ## capture effort
  eff.cap.eg <- effort %>% filter(SITE %in% 'EG')
  eff.cap.ab <- effort %>% filter(SITE %in% 'AB')
  
  if(site.keep == "El Garrapatero"){effort.hrs <- eff.cap.eg$HOUR;effort.hrs <- c(effort.hrs,as.vector(eg.net.hours))}
  if(site.keep == "Academy Bay"){  effort.hrs <- eff.cap.ab$HOUR; effort.hrs <- c(effort.hrs,as.vector(ab.net.hours))}
  names(effort.hrs) <- 2003:2018

  yr.vars <- data.frame(precipitation.year =
                          sc.num(precip.year[as.character(yr.keep)]),
                        precipitation.year.no.std =
                          (precip.year[as.character(yr.keep)]),
                        sum.rainfall =
                          sc.num(precip.btwn.samples[as.character(yr.keep)]),
                        sum.rainfall.no.std =
                          (precip.btwn.samples[as.character(yr.keep)]),
                        
                        wet.rainfall =
                          sc.num(wet.rain.our.study[wet.rain.our.study$Year.marco %in% as.character(yr.keep),"sum.prec"]),
                        wet.rainfall.no.std =
                          (wet.rain.our.study[wet.rain.our.study$Year.marco %in% as.character(yr.keep),"sum.prec"]),
                        dry.rainfall =
                          sc.num(dry.rain.our.study[dry.rain.our.study$Year.marco %in% as.character(yr.keep),"sum.prec"]),
                        dry.rainfall.no.std =
                          (dry.rain.our.study[dry.rain.our.study$Year.marco %in% as.character(yr.keep),"sum.prec"]),
                        
                        eff.cap.hour.std =
                          sc.num(effort.hrs[as.character(yr.keep)]),
                        eff.cap.hour.no.std =
                          (effort.hrs[as.character(yr.keep)]))
  
  ## *** Spline prep ***
  trait.data <- bird.data[,traits,drop=F]
  
  ## covariates
  ## cov = cbind(bird.data$PC1, # Toggle this to have only one trait
  ##             bird.data$PC2)
  
  if(is.null(K.nots)) K.nots <- max(20,min(nrow(bird.data)/4,150))

  nknotsb <- K.nots
  
  ## Toggle this to have only one trait
  R <- trait.data
  R <- R + rnorm(length(unlist(R)), mean=1e-14, sd=1e-14)
  
  knots.fields.init <- cover.design(R,nknotsb) # From fields
                                        # package, May take
                                        # some time Computes
                                        # Space-Filling
                                        # "Coverage" designs
                                        # using Swapping
                                        # Algorithm. The best
                                        # design in the form
                                        # of a matrix for the
                                        # knots minimize a
                                        # geometric
                                        # space-filling
                                        # criterion
  knots.fields <- knots.fields.init$design[,,drop=FALSE]

  Xsigma <- cbind(rep(1,nknotsb), knots.fields) # Knots matrix
  
  ## distance between the knots 
  dist.mat_all <- as.matrix(dist(knots.fields, diag=T, upper=T))

  tps.cov <- function(r) {
    f <- function(x) x^2*log(abs(x))
    r[r!=0] <- f(r[r!=0])
    r
  }
  
  ## Omega matrix [C(r) = r^2log|r|]
  OMEGA_all <- tps.cov(dist.mat_all)
  
  diffs.1_all <- outer(trait.data[,1], knots.fields[,1], '-')
  if(length(traits)==2) {
    diffs.2_all <- outer(trait.data[,2], knots.fields[,2], '-')
    dists_all <- sqrt(diffs.1_all^2 + diffs.2_all^2)
  } else {
    dists_all <- sqrt(diffs.1_all^2) # using one trait (with old code)
    # knots<-quantile(unique(trait.data[,1]), # with Crainiceanu 2005 (Bayes. A. P-Spline Reg Using BUGS)
    #                 seq(0,1,length=(nknotsb+2))[-c(1,(nknotsb+2))])
    # Z_K<-(abs(outer(trait.data[,1],knots,"-")))^3
    # OMEGA_all<-(abs(outer(knots,knots,"-")))^3
  }

  Z.1 = NULL
  Z.2 = NULL
  
  if(two.univariate.traits){
    Z_K.1 <- (abs(outer(trait.data[,1],knots.fields[,1],"-")))^3
    OMEGA_all.1 <- (abs(outer(knots.fields[,1],knots.fields[,1],"-")))^3
    ## Find the singular value decomposition of Omega 
    svd.OMEGA_all.1 <- svd(OMEGA_all.1)
    ## use SVD to obtain the matrix square root of Omega: Omega^(1/2) = U diag(sqrt(d)) V^T
    sqrt.OMEGA_all.1 <- t(svd.OMEGA_all.1$v %*% (t(svd.OMEGA_all.1$u)*sqrt(svd.OMEGA_all.1$d))) 
    ## Z = Z_kOmega^(-1/2)
    Z.1 <- t(solve(sqrt.OMEGA_all.1, t(Z_K.1) ) ) # Crainiceanu 2005 (1 trait)
    
    # Adding a second trait 
    if(length(traits)==2) {
      Z_K.2 <- (abs(outer(trait.data[,2],knots.fields[,2],"-")))^3
      OMEGA_all.2 <- (abs(outer(knots.fields[,2],knots.fields[,2],"-")))^3
      ## Find the singular value decomposition of Omega 
      svd.OMEGA_all.2 <- svd(OMEGA_all.2)
      ## use SVD to obtain the matrix square root of Omega: Omega^(1/2) = U diag(sqrt(d)) V^T
      sqrt.OMEGA_all.2 <- t(svd.OMEGA_all.2$v %*% (t(svd.OMEGA_all.2$u)*sqrt(svd.OMEGA_all.2$d))) 
      ## Z = Z_kOmega^(-1/2)
      Z.2 <- t(solve(sqrt.OMEGA_all.2, t(Z_K.2) ) ) # Crainiceanu 2005 (1 trait)
    }
  }
  ## Find the singular value decomposition of Omega 
  svd.OMEGA_all <- svd(OMEGA_all)

  ## use SVD to obtain the matrix square root of Omega: Omega^(1/2)
  ## = U diag(sqrt(d)) V^T
  sqrt.OMEGA_all <-
    t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d))) 

  ## tps.cov(dists_all) is the Z_k matrix 
  ## Z = Z_kOmega^(-1/2)
  # if(length(traits)==2) {
    Z <- t(solve(sqrt.OMEGA_all, t(tps.cov(dists_all))))
  # } else {
  #   Z <- t(solve(sqrt.OMEGA_all, t(Z_K) ) ) # Crainiceanu 2005 (1 trait)
  # }
  
  ## compute the date of first capture
  first.Xinit <- function(N, nyear, alive) {
    First <- NULL
    for (i in 1:N) {
      temp <- 1:nyear
      First <- c(First,min(temp[alive[i,]==1]))
    }
    ## inits for the states
    Xinit <- matrix(NA,nrow=N,ncol=nyear)
    for (i in 1:N) {
      for (j in 1:nyear) {
        if (j > First[i]) Xinit[i,j] <- 1
      }
    }
    return(list(First=First, Xinit=Xinit)) 
  }
  first.Xi <- first.Xinit(N=nrow(bird.data),
                          nyear=ncol(bird.data[,ch]),
                          alive=as.matrix(bird.data[,ch]))
  ## first list of inits
  init1 <- list(alive=as.matrix(first.Xi$Xinit))
  ## second list of inits
  init2 <- list(alive=as.matrix(first.Xi$Xinit))
  ## concatenate initial values
  inits.spline <- list(init1,init2)
  
  splines.vars <- list(Xinit = first.Xi$Xinit, # Xinit, spline analysis 
                       Xsigma = Xsigma, 
                       Z = Z, # Z matrix from spline analysis 
                       Z.1 = Z.1, # Z matrix from spline analysis 
                       Z.2 = Z.2, # Z matrix from spline analysis 
                       nknotsb = nknotsb, # Number of knots 
                       inits.spline = inits.spline,
                       knots.fields = knots.fields)

  
  list(X = as.matrix(bird.data[,ch]), # Needs to be a
                                      # matrix and not a
                                      # data.frame. Dataframe will
                                      # make JAGS crash.
       first = bird.data$first, # vector of first occasion
       bir.d  = bird.data,
       ind.vars = ind.vars,
       pca.recalc = pca.recalc, 
       std.mean.sd = std.mean.sd,
       yr.vars = yr.vars,
       splines.vars = splines.vars)
}

