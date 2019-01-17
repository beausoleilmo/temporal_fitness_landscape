# In this script, I calculate the logistic regressions and then extract the selection gradients. From that, I'll be able to get the values I necesitate to compute the gradients in function of the environmental variables  

# Preparation of variables and data  --------------------------------------
source("~/scripts/logit.r")
save.data = "~/Dropbox/finch_recap/v2/src/Saving RData/biotic.factors.on.survival_Andrew_meeting_changing_PCA_SCORE_for_only_fortis.RData"
# Load the appropriate files 
setwd('~/Dropbox/finch_recap/v2')
source('~/École/École - En cours!/McGill University/McGill/A - Les cours, courses/Coding (R, latex, html)/1. R/1_Functions_(handy_ones)/PCA_custom.R')
source('src/initialize.R')
# load('data/bird.data.RData', verbose=TRUE)

# Select the interesting years 
sp.list <- c('fortis')
yr = c(2004:2014,2016:2018)
jump = 1
site.list <- "El Garrapatero"
pdf.output <- FALSE
find.peaks.and.valleys = FALSE
splines <- TRUE
standard.ized = TRUE
# If you want only the linear βX not the βx and γX^2
linear.only = TRUE  
# If you want orthogonal X values: BEWARE! This is not the same as modeling an ORTHOGONAL regression 
orthogonal.x = FALSE


# prepare the variables that I'm going to record while iterating 
fit.grad.table = list()
model.list.logistic = NULL
pseudr = NULL
gofitfit = NULL
glm.model.list = NULL
yearly.number.of.id = NULL
all.ranges.x = NULL
list.min.fit = NULL
list.min.fit.trait = NULL
original.x = NULL
midd.list = NULL
my.eco.evo.df = NULL
full.data = NULL
final.df.new.analysis = NULL
newlist = list()
oldxlist = list()
oldzlist = list()
old_beak_L_list = list()
old_beak_W_list = list()
old_beak_D_list = list()
old_band_list = list()
model.list = list()
survived.list = list()

# These were the relatively good values of lambda used 
# exp.lambda = exp(seq(-8,-1,by =2))
exp.lambda = exp(mean(c(-4,-5,-4,-4,-4,-4,-4,-4,-3,-6,-13,0)))
# exp.lambda = exp(-2)
load('data/bird.data.RData', verbose=TRUE)
dput(names(bird.data))
autopairs(x = bird.data[,c("Site", "Sex0","Site", "Sex0","Year",  "Species1",
                           "Tarsus", "Wing.Chord", "Mass",  
                           "MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth",
                           "PC1", "PC2" #"pc1.score.per.sp", "pc2.score.per.sp", "pc1.score.per.sp.eg", "pc2.score.per.sp.eg", "PC.body1", "PC.body2"
)], main = "Traits all ground finches at EG and AB")
autopairs(x = bird.data[bird.data$Species1 == "fortis", # & bird.data$Site == "El Garrapatero",
                        c("Site", "Sex0","Site", "Sex0","Year",  "Species1",
                          "Tarsus", "Wing.Chord", "Mass",
                          "MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth",
                          "PC1", "PC2" #"pc1.score.per.sp", "pc2.score.per.sp", "pc1.score.per.sp.eg", "pc2.score.per.sp.eg", "PC.body1", "PC.body2"
                        )], main = "Traits G. fortis at EG and AB")
cor(bird.data[bird.data$Species1=="fortis","Mass"], 
    bird.data[bird.data$Species1=="fortis","MedianBeakDepth"])
cor(bird.data[bird.data$Species1=="fortis","Mass"], 
    bird.data[bird.data$Species1=="fortis","MedianBeakLength"])
cor(bird.data[bird.data$Species1=="fortis","Mass"], 
    bird.data[bird.data$Species1=="fortis","MedianBeakWidth"])
cor(bird.data[bird.data$Species1=="fortis","Mass"], 
    bird.data[bird.data$Species1=="fortis","PC1"])
plot((bird.data[bird.data$Species1=="fortis","Mass"])^(1/3), 
     bird.data[bird.data$Species1=="fortis","PC1"])
cor(bird.data$Mass,bird.data$MedianBeakDepth)
cor(bird.data$Mass,bird.data$MedianBeakLength)
cor(bird.data$Mass,bird.data$MedianBeakWidth)
cor(bird.data$Mass,bird.data$Tarsus)
# bird.fort = prep.data(sp.keep = sp.list,
#                       yr.keep = c(2003:2014,2015,2016:2018),  
#                       site.keep = site.list,
#                       traits = "PC1")
if(pdf.output){
  pdf("~/Desktop/my.pca.just.fortis.pdf",height = 6,width = 15)
  par(mfrow=c(1,3))
  res.pca.for2 = vegan::rda(bird.data[bird.data$Species1=="fortis", 
                                      c("MedianBeakLength", "MedianBeakWidth", 
                                        "MedianBeakDepth")])
  res.pca.for3 = vegan::rda(bird.data[bird.data$Species1=="fortis" & bird.data$Site=="El Garrapatero",
                                      c("MedianBeakLength", "MedianBeakWidth", 
                                        "MedianBeakDepth")])
  nrow(bird.data[bird.data$Species1=="fortis" & bird.data$Site=="El Garrapatero",
                 c("MedianBeakLength", "MedianBeakWidth", 
                   "MedianBeakDepth")])
  bd.EG =  bird.data[bird.data$Site=="El Garrapatero",]
  res.pca.all.sp.EG= vegan::rda(bd.EG[,c("MedianBeakLength", 
                                         "MedianBeakWidth", 
                                         "MedianBeakDepth")])
  bd.EG =  bird.data[bird.data$Site=="El Garrapatero",]
  res.pca.all.sp.EG= vegan::rda(bd.EG[,c("MedianBeakLength", 
                                         "MedianBeakWidth", 
                                         "MedianBeakDepth")])
  
  length(bird.data[bird.data$Species1=="fortis","PC1"])
  table(bird.data$Species1,bird.data$Year)
  tab.sp.eg= table(bird.data[bird.data$Site=="El Garrapatero","Species1"], bird.data[bird.data$Site=="El Garrapatero","Year"])
  write.csv(tab.sp.eg,"~/Desktop/nb.sp.eg.csv")
  table(bird.data[bird.data$Species1=="fortis",]$Year)
  table(bird.data[bird.data$Species1=="fortis"& bird.data$Site=="El Garrapatero",]$Year)
  table(bird.data[bird.data$Species1=="fortis"& bird.data$Site=="Academy Bay",]$Year)
  length(bird.data[bird.data$Species1=="fortis" & bird.data$Site=="El Garrapatero","PC1"])
  ores = mixtools::normalmixEM(bird.data[bird.data$Species1=="fortis","PC1"],
                               # mu = c(center2,center1),
                               sigma = NULL, 
                               mean.constr = NULL, sd.constr = NULL,
                               epsilon = 1e-15, maxit = 1000, maxrestarts=50, 
                               # verb = TRUE, 
                               fast=FALSE, ECM = FALSE,
                               arbmean = TRUE, arbvar = TRUE)
  
  par(mfrow=c(1,2))
  res.pca.all.sp.EG$CA$u[,1] = -res.pca.all.sp.EG$CA$u[,1]
  res.pca.all.sp.EG$CA$v[,1] = -res.pca.all.sp.EG$CA$v[,1]
  custom_pca(res.pca.all.sp.EG, centered = TRUE, 
             vector.to.colour.col = bd.EG$Species1,
             vector.to.colour.pch = bd.EG$Species1,
             col.label = "black",
             shape.points = c(19,22,23,24),
             label.the.arrows = TRUE)
  summary(res.pca.all.sp.EG)
  round(scores(res.pca.all.sp.EG, choices = 1:3, display = "species", scaling = 0),2)
  scores(res.pca.all.sp.EG, choices = 1:3, display = "species")
  
  res.pca.for$CA$u[,2] = -res.pca.for$CA$u[,2]
  # res.pca.for$CA$u[,1] = -res.pca.for$CA$u[,1]
  # res.pca.for$CA$v[,2] = -res.pca.for$CA$v[,2]
  # res.pca.for$CA$v[,1] = -res.pca.for$CA$v[,1]
  custom_pca(res.pca.for, centered = TRUE,
             ordiellipse = TRUE,shape.points = 21,
             ordiellipse.vector = ores$posterior[,1]>.95,
             col.label = "black",
             label.the.arrows = FALSE,
             vector.names = c("Median Beak Length",
                              "Median Beak Width",
                              "Median Beak Depth"))
  
  
  
  custom_pca(res.pca.for2, centered = TRUE,
             ordiellipse = TRUE,shape.points = 21,
             # ordiellipse.vector = ores$posterior[,1]>.95,
             col.label = "black",
             label.the.arrows = FALSE,
             vector.names = c("Median Beak Length",
                              "Median Beak Width",
                              "Median Beak Depth"))
  
  pc1.raw=scores(res.pca.for3,choices = 1)
  oorees = mixtools::normalmixEM(pc1.raw$sites,
                                 # mu = c(center2,center1),
                                 sigma = NULL, 
                                 mean.constr = NULL, sd.constr = NULL,
                                 epsilon = 1e-15, maxit = 1000, maxrestarts=50, 
                                 # verb = TRUE, 
                                 fast=FALSE, ECM = FALSE,
                                 arbmean = TRUE, arbvar = TRUE)
  res.pca.for3$CA$u[,2] = -res.pca.for3$CA$u[,2]
  res.pca.for3$CA$v[,2] = -res.pca.for3$CA$v[,2]
  custom_pca(res.pca.for3, centered = TRUE,
             ordiellipse = TRUE,shape.points = 21,
             ordiellipse.vector = oorees$posterior[,1]>.95,
             col.label = "black",
             label.the.arrows = FALSE,
             vector.names = c("Median Beak Length",
                              "Median Beak Width",
                              "Median Beak Depth"))
  summary(res.pca.for3)
  round(scores(res.pca.for3, choices = 1:3, display = "species", scaling = 0),2)
  scores(res.pca.for3, choices = 1:3, display = "species")
  
  #                    PC1     PC2    PC3
  # MedianBeakLength 4.223  2.0071  0.166
  # MedianBeakWidth  4.637 -0.4482 -1.391
  # MedianBeakDepth  6.200 -1.0319  0.927
  
  
  # res.pcabs is found in load('data/bird.data.RData', verbose=TRUE)
  res.pcabs$CA$u[,1] = -res.pcabs$CA$u[,1]
  custom_pca(res.pcabs, centered = TRUE, 
             vector.to.colour.col = bird.data$Species1,
             vector.to.colour.pch = bird.data$Species1,
             col.label = "black",
             shape.points = c(19,22,23,24),
             label.the.arrows = TRUE)
  
  summary(res.pcabs)
  summary(res.pca.all.sp.EG)
  dev.off()
}


if(pdf.output){
  pdf(paste0("~/Desktop/my.fit.land.short2_jump",
             jump,"_",gsub('([[:punct:]])|\\s+','_',site.list),".pdf"),
      height = 7,
      width = 8)
}

# Function finding maximum and minimum of the fitness function  -----------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
par(mfrow=c(1,1))
# For loop that will calculate the GAM, show the landscape and will let you select what is the maximum and minimum of the function 
if(find.peaks.and.valleys){
  for(i in 1:c(length(yr)-1)){
    yr.list <- c(yr[i],yr[i+jump])
    
    ## If you want to try only fuliginosa, use only this species in the sp.list 
    mdat <- prep.data(sp.keep = sp.list,
                      center = FALSE,
                      yr.keep = yr.list,
                      site.keep = site.list,
                      flip.pc1 = FALSE,
                      adults.only = FALSE, # If true, it'll keep only females and males (remove juveniles) 
                      gr.sel = NULL,# "small" or "big", this is to create 2 groups in the fortis population. That way, I could see if selection is acting differently on the 2 fortis groups, # applicable for fortis only. This will let you select "big" or "small" morph
                      valley.fortis = FALSE, # If this is true, this term will find the valley (between the small morph of fortis and the big morph of fortis). You can run a model on these 2 after that
                      valley.fortis.fuli = FALSE, # If this is true, this term will find the valley (between fuliginosa and the small morph of fortis). You can run a model on these 2 after that
                      peaks.fortis.fuli = FALSE, # If this is true, this term will find the 2 peaks (of fuliginosa and the one for the small morph of fortis). You can run a model on these 2 after that
                      two.univariate.traits = FALSE, # For spline, # For spline only. This is computing PC1 and PC2 as 2 independent univariate traits.
                      keep.resid = FALSE, # FALSE if want to keep everything, # If this is true, it's going to keep the individuals that are putatively resident (==1)
                      splines = splines,
                      keep.last.year.data = FALSE,
                      gam.analysis = TRUE,
                      recalculate.pca = FALSE,
                      pca.per.sp.only = TRUE,
                      K.nots = 5, # Manually select the number of knots  
                      traits = c('PC1'))#, 'PC2')) # This is only for the spline 
    if(!is.null(mdat$pca.recalc)){
      custom_pca(mdat$pca.recalc, centered = TRUE)
    }
    custom_pca(res.pca.for, centered = TRUE)
    # This is to calcualte the GAM 
    library(mgcv)
    
    # getting response variable and the explanatory variable 
    y = as.vector(mdat$X[,2])
    x = c(mdat$ind.vars$pc1) # Works for everything execpt 2007-2008
    mbd = c(mdat$ind.vars$mbd) 
    mbl = c(mdat$ind.vars$mbl) 
    mbw = c(mdat$ind.vars$mbw) 
    band = as.character(mdat$ind.vars$band) 
    year.var = rep(yr.list[2],length(mdat$ind.vars$pc1))
    # x = c(mdat$ind.vars$mbd) # Works for most years 
    # x = c(mdat$ind.vars$mbw) # Works only for 2004-2005 and maybe 2008-2009
    # x = c(mdat$ind.vars$mbl) # Works only for 2005-2006 and maybe 2006-2007, perhaps 2008-2009, 2009-2010
    # x = c(mdat$ind.vars$mass) # Effectively works for only 2004-2005 and probably 2008-2009
    mydata = data.frame(x,y)
    dat.for.comparison.analysis = data.frame(x,y,year.var)
    full.data = c(full.data,list(dat.for.comparison.analysis))
    # Fitting the GAM 
    z <- gam(y ~ s(x), 
             data = mydata, 
             family = binomial(link = "logit"), # Logistic link 
             sp = exp.lambda,
             method = "GCV.Cp")
    
    model.list = c(model.list,list(z))
    
    # Main model 
    # plot(z, se = 1, seWithMean = TRUE, rug = FALSE, shift = mean(predict(z)),
    #      ylim = c(0,1),
    #      trans = function(x){exp(x)/(1+exp(x))}, 
    #      main = paste("GAM ±2SE binom","yr", i+2003, i+2003+jump, sep = " "))  # binomial data
    yr1 = substr(i+2003,3,4)
    yr2 = substr(i+2003+jump,3,4)
    # Same thing but with all points 
    # lambb.z=round(log(z$sp))
    
    # Show Gam with only points of individuals present 
    # plot(z$fitted.values~x, 
    #      cex = .8, pch = 21, 
    #      bg = "black", col = "black")
    # 
    oldxlist = c(oldxlist, list(x))
    oldzlist = c(oldzlist, list(z$fitted.values))
    old_beak_L_list = c(old_beak_L_list, list(mbd))
    old_beak_W_list = c(old_beak_W_list, list(mbl))
    old_beak_D_list = c(old_beak_D_list, list(mbw))
    old_band_list = c(old_band_list, list(band))
    
    survived.list = c(survived.list, list(y))
    
    lambb.z = #lamb=
      round(log(exp.lambda))
    
    # Getting new x that is spaced evenly respecting the GAM function (fitness function) 
    newx <- seq(from = min(mydata$x), 
                to = max(mydata$x), 
                length.out = 2000)
    # Using the model to generate the new response variable 
    z1 <- predict(z, 
                  newdata=list(x = newx), 
                  se.fit = TRUE)
    # Function that will transforme the Y values from linear scale to [0,1]
    invlogit <- function(x){exp(x)/(exp(x) + 1)}
    # This is the actual transformed data 
    yhat <- invlogit(z1$fit)
    upper <- invlogit(z1$fit + z1$se.fit)
    lower <- invlogit(z1$fit - z1$se.fit)
    
    # satisfied...= "n"
    # while(satisfied... != "y"){
    # Here is the plot of the fitness function using the evenly spaced data (newx) and the response to it (yhat)
    plot(newx, yhat, type="l", 
         ylim = c(0,1), 
         xlab = "PC1",
         main = bquote(atop("GAM±1SE",
                            lambda * " = "*.(lambb.z) * 
                              ",  yr = "*.(yr1) *"-"* .(yr2))))
    # Adding error 
    lines(newx, upper, lty = 2)
    lines(newx, lower, lty = 2)
    # par(new = TRUE)
    # dtrait = density(x)
    # plot(dtrait, xlim = range(x), col = "blue")
    # These values are going to be emplemented in order to find the maximum (2 values) and minimum (1 value) of the fitness function 
    
    # find the minimum of the fitness function from the gam by clicking on BOTH sides of the highest visible peak and  the minimum value between the 2 peaks (valley) in the GAM  
    midd  = locator(n = 2)
    # Starting from the left side, click right around the maximum of the fitness function. This will find the maximum value 
    peak  = locator(n = 4)
    # midd=NULL
    # This will extract the X values for the selections 
    midd = midd$x
    peak = peak$x
    # midd = readline(prompt="Enter middle to find local minima: ")
    
    # Here is the actual function that will find the maximum and minium 
    local.min  = min(yhat[newx > midd[1] & newx < midd[2]])
    local.max  = max(yhat[newx > midd[1] & newx < midd[2]])
    
    # make a database usign the 2 peaks 
    local.max.peak1  = newx[which(yhat==max(yhat[newx > peak[1] & newx < peak[2]]))]
    local.max.peak2  = newx[which(yhat==max(yhat[newx > peak[3] & newx < peak[4]]))]
    
    # Draw a diagnostic line to see that you've selected the mximium of an obserbed peak and the minimum of the valley 
    abline(h = c(local.min, local.max), lty = 3)
    # This is showing the 2 maximum values of the 2 peaks 
    abline(v = c(local.max.peak1, local.max.peak2), lty = 3, col = "red")
    
    # This si a metric of the fitness function itself. That means that the highest fitness- the lowest fitness would be the fitness differential (if it exists)
    midd.list = c(local.max - local.min)
    
    # Make a record of all the expected response varaible from the evenly spaced X 
    newlist <- c(newlist, list(yhat))
    
    # }
    
    # Make a database for all these new varaibles 
    my.eco.evo.df=rbind(my.eco.evo.df,data.frame(mid = midd.list,
                                                 # yr = paste(yr.list,collapse = '_'), 
                                                 yr1 = yr.list[1],
                                                 yr2 = yr.list[2], 
                                                 local.max.peak1 = local.max.peak1,
                                                 local.max.peak2 = local.max.peak2,
                                                 preci.yr1 = mdat$yr.vars$precipitation.year.no.std[1], 
                                                 preci.yr2 = mdat$yr.vars$precipitation.year.no.std[2],
                                                 sum.preci.yr1 = mdat$yr.vars$sum.rainfall.no.std[1], 
                                                 sum.preci.yr2 = mdat$yr.vars$sum.rainfall.no.std[2]))
    
  }# End of for(i in 1:c(length(yr)-1)){
  
  # dev.off()
  
  eff = structure(c(36, 140, 212, 120, 52, 56, 132, 300, 128, 120, 128, 
                    104, 30.2565802161513, 81.4673085182818, 81.4673085182818, 117.180771285717
  ), .Names = c("2003", "2004", "2005", "2006", "2007", "2008", 
                "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", 
                "2017", "2018"))
  my.eff = eff[names(eff)%in% yr]
  
  par(mfrow=c(2,2))
  heights = my.eco.evo.df[,1]
  my.eco.evo.df$good = c(1,1,0,1,1,1,0,0,0,0,0,0,0)
  my.eco.evo.df$good = c(1,1,1,1,1,1,1,0,0,0,0,0,0)
  my.eco.evo.df$good = rep(1, nrow(my.eco.evo.df))
  my.eco.evo.df$peak.height = as.factor(c("l","r","r","r",
                                          "l","r","l","r",
                                          "l","m","m","m",
                                          "m"))
  my.eco.evo.df$two.peaks = c(0,1,1,1,
                              1,1,1,0,
                              0,0,0,0,
                              0)
  my.eco.evo.df[my.eco.evo.df$peak.height == "m",""]
  # plot(1:3, col=1:3)
  # Keep only the relevant year combination
  my.eco.evo.df$effort = my.eff[1:c(length(my.eff)-1)]
  my.eco.evo.df$effort = my.eco.evo.df$effort/max(my.eco.evo.df$effort)
  my.df = my.eco.evo.df[my.eco.evo.df$good == 1,]
  
  fct = (mid)~(sum.preci.yr2)
  lm.out1 = lm(fct,data = my.df)
  summary(lm.out1)
  plot(fct,data = my.df, 
       cex = my.df$effort, 
       pch =21, 
       bg = my.df$peak.height, 
       col = my.df$peak.height,
       ylab = "Valley depth (in fitness difference)",
       xlab = "Cumulative precipitation");abline(lm.out1)
  text(x = my.df$sum.preci.yr1,
       y = my.df$mid,
       labels = my.df$yr2, cex = .6, pos = 3)
  # lines(x = (my.df$sum.preci.yr2),
  #       y = (my.df$mid), lwd = .1)
  
  fct = log((mid))~log((sum.preci.yr2))
  lm.out2=lm(fct,data = my.df)
  summary(lm.out2)
  plot(fct, data = my.df, 
       cex = my.df$effort, 
       pch =21, 
       bg = my.df$peak.height, 
       col = my.df$peak.height, 
       ylab = "ln(Valley depth) (in fitness difference)",
       xlab = "ln(Cumulative precipitation)");abline(lm.out2)
  text(x = log(my.df$sum.preci.yr1),
       y = log(my.df$mid),
       labels = my.df$yr2, cex = .6, pos = 3)
  # lines(x = log(my.df$sum.preci.yr2),
  #       y = log(my.df$mid), lwd = .1)
  # plot(lm.out2)
} # End of if(find.peaks.and.valleys){


# Save finding peaks ------------------------------------------------------
# Save.once
if(find.peaks.and.valleys){
  save(my.eco.evo.df,
       lm.out1,
       lm.out2,
       newlist,
       newx,
       model.list,
       oldxlist,
       oldzlist,
       old_beak_L_list,
       old_beak_W_list,
       old_beak_D_list,
       old_band_list,
       survived.list,
       full.data,
       file = save.data)
}
# dev.off()


# Load the data  ----------------------------------------------------------
load(save.data, 
     verbose = TRUE)

# Here I wanted to test if it is possible to find differences between the variance of 2 distribution of beaks in 2 different years for the same beak type (for example, compare the beak distribution variance in 2005 for the big beak to the variance in distribution of the beak beak in 2006). => I compared 2005-2006 for beak beak and the variance in beak measurement changes! 
ores.beak = list()
for (i in 1:length(full.data)) {
  ores.beak.tmp = mixtools::normalmixEM(full.data[[i]]$x,
                                        # mu = c(center2,center1),
                                        sigma = NULL, 
                                        mean.constr = NULL, sd.constr = NULL,
                                        epsilon = 1e-15, maxit = 1000, maxrestarts=50, 
                                        # verb = TRUE, 
                                        fast=FALSE, ECM = FALSE,
                                        arbmean = TRUE, arbvar = TRUE)
  ores.beak = c(ores.beak,list(ores.beak.tmp))
}

find.max = apply(ores.beak[[1]]$posterior,1,which.max)
mixturemodel.out.groups1 = data.frame(pc1 = ores.beak[[1]]$x,
                                      gr = find.max)
find.max = apply(ores.beak[[2]]$posterior,1,which.max)
mixturemodel.out.groups2 = data.frame(pc1 = ores.beak[[2]]$x,
                                      gr = find.max)
compar= data.frame(pc1 = c(mixturemodel.out.groups1[mixturemodel.out.groups1$gr==1,"pc1"],
                           mixturemodel.out.groups2[mixturemodel.out.groups2$gr==1,"pc1"]),
                   gr = c(rep(x = "1",length(which(mixturemodel.out.groups1$gr ==1))),
                          rep(x = "2",length(which(mixturemodel.out.groups2$gr ==1)))))

# http://www.sthda.com/english/wiki/compare-multiple-sample-variances-in-r#compute-levenes-test-in-r
# Statistical tests for comparing variances
# There are many solutions to test for the equality (homogeneity) of variance across groups, including:
#   F-test: 
# Compare the variances of two samples. The data must be normally distributed.
#   Bartlett’s test: 
# Compare the variances of k samples, where k can be more than two samples.  The data must be normally distributed. The Levene test is an alternative to the Bartlett test that is less sensitive to departures from normality.
#   Levene’s test: 
# Compare the variances of k samples, where k can be more than two samples. It’s an alternative to the Bartlett’s test that is less sensitive to departures from normality.
#   Fligner-Killeen test: 
# a non-parametric test which is very robust against departures from normality.


bartlett.test(pc1 ~ gr,data=compar)
# If not normally distributed 
car::leveneTest(pc1 ~ gr,data=compar)

# plot(mixturemodel.out.groups$pc1, col = mixturemodel.out.groups$gr)
par(mfrow=c(3,4))
for (i in 1:12) {
  plot(ores.beak[[i]],whichplots = 2, main2 = c(2005:2018)[i])
}
par(mfrow=c(1,1))


my.eco.evo.df

par(mfrow= c(2,1))
mydf = data.frame(y = newlist[[1]], x = newx)
where.peak1 = which(mydf$y == max(mydf$y[mydf$x>my.eco.evo.df$local.max.peak1[1]-.01  & 
                                           mydf$x<my.eco.evo.df$local.max.peak1[1]+.1]))
# where.peak1 = 1

# newlist[[1]] > my.eco.evo.df$local.max.peak1[1]
# my.eco.evo.df$
# round(newx,4) >= round(my.eco.evo.df$local.max.peak1[1],4)
sub.newx = newx[where.peak1:nrow(mydf)]
x = sub.newx
q.out = lm(mydf$y[where.peak1:nrow(mydf)]~poly(x,2))
summary(q.out)

# with old data -----------------------------------------------------------
# This is an attempt at looking at the QUADRATIC curves. But this is mathematically wrong and not advised. Andrew hasn't done that 
# plot(mydf$y~newx, ylim=c(0,1), xlim=range(newx), type ="l")
# curve(predict(q.out,data.frame(x=x)), col = "red",lwd=2,add=TRUE)
# plot(mydf$y[where.peak1:nrow(mydf)]~newx[where.peak1:nrow(mydf)], ylim=c(0,1), xlim=range(newx), type ="l")
# curve(predict(q.out,data.frame(x=x)), col = "red",lwd=2,add=TRUE)
#
# mydf = data.frame(y = oldzlist[[1]], x = oldxlist[[1]], survived.list[[1]])
# mydf=mydf[order(mydf$x),]
# mydf = mydf[mydf$survived.list..1.. == 1,]
# where.peak1=which(mydf$y ==max(mydf$y[mydf$x>my.eco.evo.df$local.max.peak1[1]-.03  & 
#                                         mydf$x<my.eco.evo.df$local.max.peak1[1]+.03]))
# # where.peak1 = 2
# x = mydf$x[where.peak1:nrow(mydf)]
# q.out = lm(mydf$y[where.peak1:nrow(mydf)]~poly(x,2))
# 
# plot(mydf$y~mydf$x, 
#      ylim=c(0,1), 
#      xlim=range(newx), type ="l", ylab= "Survival", xlab ="PC1", 
#      main ="fitting quadratic function \n to interesting section")
# curve(predict(q.out,data.frame(x=x)), col = "red",lwd=2,add=TRUE)
# plot(mydf$y[where.peak1:nrow(mydf)]~x,
#      ylim=c(0,1), xlim=range(mydf$x), type ="l",
#      ylab= "Survival", xlab ="PC1")
# points(mydf$y[where.peak1:nrow(mydf)]~x, cex = .7, pch = 21, bg = "black")
# curve(predict(q.out,data.frame(x=x)), col = "red",lwd=2,add=TRUE)
# summary(q.out)


# These are the points that I need! 
par(mfrow=c(2,3))
sample.size = NULL
my.eco.evo.df = my.eco.evo.df[1:7,]


for(i in 1:nrow(my.eco.evo.df)){
  ah.analysis = data.frame(y =  oldzlist[[i]], 
                           x =  oldxlist[[i]], 
                           mbl = old_beak_L_list[[i]],
                           mbw = old_beak_W_list[[i]],
                           mbd = old_beak_D_list[[i]],
                           band = old_band_list[[i]],
                           app.surv = survived.list[[i]])
  ah.analysis = ah.analysis[order(ah.analysis$x),]
  where.peak1 = which(ah.analysis$y == max(ah.analysis$y[ah.analysis$x>my.eco.evo.df$local.max.peak1[i]-.11  & 
                                                           ah.analysis$x<my.eco.evo.df$local.max.peak1[i]+.11]))
  where.peak2 = which(ah.analysis$y == max(ah.analysis$y[ah.analysis$x>my.eco.evo.df$local.max.peak2[i]-.11  & 
                                                           ah.analysis$x<my.eco.evo.df$local.max.peak2[i]+.11]))
  
  # ah.analysis = ah.analysis[!(ah.analysis$x >.8 & ah.analysis$x <1 &ah.analysis$app.surv == 1),]
  sub = ah.analysis[where.peak1:where.peak2,]
  # plot(sub)
  
  original.x = c(original.x,list(sub$x))
  my.x = sub$x
  my.mbl = sub$mbl
  my.mbw = sub$mbw
  my.mbd = sub$mbd
  my.band = sub$band
  y = sub$app.surv
  # plot(y~x)
  yea.var=rep(my.eco.evo.df[i,"yr2"],length(y))
  new.analysis.all.years = data.frame(my.x,y,yea.var,my.mbl,my.mbw,my.mbd,my.band)
  final.df.new.analysis = c(final.df.new.analysis,list(new.analysis.all.years))
  # Preparing the data X (quadratic or not, orthongonal or not)
  
  stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
  if (standard.ized) {
    x.se = stderr(my.x)
    my.x = scale(my.x)
    x.mean = attr(my.x,"scaled:center")
    x.sd = attr(my.x,"scaled:scale")
    # Retransformed data 
    x.st = data.frame(x = my.x*x.sd+x.mean, 
                      x2 = (my.x*x.sd+x.mean)^2)
  }
  
  if (linear.only) {
    x = data.frame(x = my.x)
  } else if(orthogonal.x) {
    x = data.frame(x =  poly(my.x,degree = 2,raw = FALSE)[,1], 
                   x2 = poly(my.x,degree = 2,raw = FALSE)[,2])  
  } else {
    x = data.frame(x = my.x, 
                   x2 = my.x^2)
  }
  
  all.ranges.x = c(all.ranges.x,list(range(x$x*x.sd+x.mean)))
  
  
  # According to endler (p.255), if the covariation between X and X^2 is 0, there is no skweness and therefore the univariate regression and multivariate regression coefficients are the same 
  library(moments)
  kurtosis(x$x)
  skewness(x$x) # 
  cov(x$x,(x$x-mean(x$x))^2)
  # When standardized 
  # 0.08425317
  
  galtonskew.proc <- function(x){
    #
    #  Compute Galton's skewness measure for x
    #  NOTE: this procedure assumes no x values are missing
    #
    quarts <- as.numeric(quantile(x, probs = c(0.25, 0.5, 0.75)))
    num <- quarts[1] + quarts[3] - 2*quarts[2]
    denom <- quarts[3] - quarts[1]
    gskew <- num/denom
    gskew
  }
  galtonskew.proc(x$x)
  
  sample.size = c(sample.size,
                  list(table(y)))
  library("R2ucare")
  ch = cbind(x= rep(1,length(y)),y)
  # table(apply(ch, 1, function(x) {paste(x,collapse = "")}))
  freq = rep(1,length(y))
  # overall_CJS(ch,freq)
  # test2cl(X = ch,freq)
  # test2ct(X = ch,freq = freq)
  # test3Gsm(ch, freq)
  # test3Gsr(ch, freq)
  # test3Gwbwa(ch, freq)
  # test3sm(ch, freq)
  # test3sr(ch, freq)
  
  model = glm(y ~ ., 
              family=binomial, 
              data = x)
  glm.model.list = c(glm.model.list,
                     list(model))
  
  if (standard.ized) {
    # This uses the original data since x.st is the retransformation from the scaled to the original data 
    model.st = glm(y ~ ., 
                   family=binomial, 
                   data = x.st)
  }
  model
  pseudoR2 <- (model$null.deviance - model$deviance) / model$null.deviance
  # Same calculation McFadden’s Pseudo-R^2
  # pR2 = 1 - model$deviance / model$null.deviance # works for glm
  
  print(pseudoR2)
  # PseudoR2(model)
  library("vcdExtra")
  gofit = HLtest(model)
  cat("The p-value of the GOF is:",round(gofit$p.value,2),"\n")
  
  pseudr = c(pseudr,pseudoR2)
  gofitfit = c(gofitfit,list(gofit))
  
  
  model.list.logistic = rbind(model.list.logistic,
                              as.table(coefficients(model)))
  sumglm = summary(model)
  
  
  # Compute Gradients -----------------------------------------
  #### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
  # http://www.public.iastate.edu/~fjanzen/software/regress.htm
  
  # To get better estimate of the linear selection gradient: run the models separately  (w = α + βx and then w = α + βx + γx^2)
  if (length(dim(x)) != 2) dim(x) <- c(length(x),1)
  # For the linear comparison 
  alin <- lsfit(x,y) # Find the Least Squares Fit
  blin <- ls.diag(alin) # Compute Diagnostics for lsfit Regression Results
  # For the logistic  comparison 
  algt <- model #glm(y ~ ., family = binomial(link = "logit"), data=x) # logistic regression
  blgt <- summary(algt)
  
  n <- nrow(x)
  np <- ncol(x)
  rw <- n/sum(y) # Relative fitness 
  mn <- apply(x,2,mean) # mean of each trait 
  md <- apply(x,2,median) # median of each trait 
  sd <- sqrt(apply(x,2,var))[1] # standard deviation, same as apply(x,2,sd) 
  mean.z <- apply(x,2,mean)[1] # mean trait
  
  # Why estimate linear and full (linear, quadratic, and correlational) selection coefficients separately?
  # https://biology.stackexchange.com/questions/24728/why-estimate-linear-and-full-linear-quadratic-and-correlational-selection-co
  
  # About orthogonal regressions
  # orthogonal regression. This is where we fit a regression line so that we minimize the sum of the squares of the orthogonal (rather than vertical) distances from the data points to the regression line.
  
  # https://stats.stackexchange.com/questions/13152/how-to-perform-orthogonal-regression-total-least-squares-via-pca
  v <- prcomp(cbind(x,y))$rotation
  beta <- v[2,1]/v[1,1]
  # See also 
  # library(onls)
  # For the linear comparison 
  
  c1 <- alin$coef[2:(np+1)]   *sd*rw # this is rescaling the coefficients with sd of the traits AND relative fitness; (2:(np+1) is removing the intercept)
  c2 <- blin$std.err[2:(np+1)]*sd*rw # this is rescaling the std.err with sd of the traits AND w 
  c1[2] = c1[2]*sd
  # c1[2] = c1[2]*mean.z
  c2[2] = c2[2]*sd
  # c2[2] = c2[2]*mean.z
  c3 <- 2*(1-pt(q = abs(c1)/c2,
                df = n-np-1)) # pt: The Student t Distribution; q is vector of quantiles; n-np-1 are the degrees of freedom; 2* is because it's 2 tailed 
  # For the logistic  comparison (alpha = coef1*sdev1;)
  c4 <- blgt$coefficients[2:(np+1),1]*sd # Estimates of the GLM (logistic regression), scaled by the sd 
  c4.1 <- blgt$coefficients[2:(np+1),1]*mean.z # Estimates of the GLM (logistic regression), scaled by the sd 
  c4[2] = c4[2]*sd
  c4.1[2] = c4.1[2]*mean.z
  # SE (se = serr*sdev1;)
  c5 <- blgt$coefficients[2:(np+1),2]*sd # SE of the GLM (logistic regression), scaled by the sd 
  c5.1 <- blgt$coefficients[2:(np+1),2]*mean.z # SE of the GLM (logistic regression), scaled by the sd 
  c5[2] = c5[2]*sd
  c5.1[2] = c5.1[2]*mean.z
  
  c6 <- 2*(1-pnorm(q = abs(c4)/c5, 
                   mean = 0, 
                   sd = 1, 
                   lower.tail = TRUE, 
                   log.p = FALSE)) # pnorm: Normal Distribution to find the 2 tailed P-value (2* and abs() function)
  c6.1 <- 2*(1-pnorm(q = abs(c4.1)/c5.1, 
                     mean = 0, 
                     sd = 1, 
                     lower.tail = TRUE, 
                     log.p = FALSE)) # pnorm: Normal Distribution to find the 2 tailed P-value (2* and abs() function)
  
  # Get the predicted values on the LINEAR SCALE; verify like this range(z)
  z <- blgt$coefficients[1,1] + # intercept 
    as.matrix(x) %*% blgt$coefficients[2:(np+1),1] # α_0 + α_1*x_1 + α_i x_i + ... + α_n + x_n
  # Transform 
  invlogit <- function(x) {exp(x)/(1+exp(x))}
  pr <- invlogit(z) # transormation to the logistic scale [0,1]; verify like this range(pr)
  # (bavgrad = coef1*sdev1*ag/rf;) => bavgrad = c4 * mean(pr*(1-pr)) * rw
  c7 <- mean(pr*(1-pr))*c4*rw # pr*(1-pr) refers to {\hat{p}*(1-\hat{p})}, c4 contains the estimates from the GLM and rw is the relative fitness. This line is an approximation of the average selection gradient
  c7.1 <- mean(pr*(1-pr))*c4.1*rw # pr*(1-pr) refers to {\hat{p}*(1-\hat{p})}, c4 contains the estimates from the GLM and rw is the relative fitness. This line is an approximation of the average selection gradient
  c8 <- mean(pr*(1-pr))*c5*rw 
  c8.1 <- mean(pr*(1-pr))*c5.1*rw 
  
  # Multiply by 2 the estimates 
  c1[2] = c1[2]*2
  c2[2] = c2[2]*2
  c4[2] = c4[2]*2
  c5[2] = c5[2]*2
  c7[2] = c7[2]*2
  c7.1[2] = c7.1[2]*2
  c8[2] = c8[2]*2
  c8.1[2] = c8.1[2]*2
  
  ag = mean(pr*(1-pr))
  cag = rep(ag, length(c1))
  x.mean = rep(x.mean, length(c1))
  x.sd.rep = rep(x.sd, length(c1))
  x.se.rep = rep(x.se, length(c1))
  
  Relfit = rep(1/rw, length(c1))
  
  c0 = blgt$coefficients[2:(np+1),1]
  
  # Create the table 
  tab <- cbind(c0, c1, c2, c3,
               c4, c5, c6,
               c7, c8, 
               cag, Relfit, sd, mean.z,
               c4.1, c5.1, c6.1,
               c7.1, c8.1,x.mean,
               x.sd,x.se)
  dimnames(tab)[[2]] <- c("coeff", "lin.beta", "lin.se", "lin.p",
                          "lgt.alph", "lgt.se", "lgt.p",
                          "bavggrd", "bavggrds", 
                          "cag", "Relfit", "sd", "mean",
                          "lgt.alph.mean.st", "lgt.se.mean.st", "lgt.p.mean.st",
                          "bavggrd.mean.st", "bavggrds.mean.st","x.mean","x.sd","x.se"
  )
  if (linear.only) {
    tab = tab[seq(1,nrow(tab),by = 2),,drop=FALSE]
  }
  fit.grad.table = c(fit.grad.table,list(as.data.frame(tab)))
  #### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
  
  
  
  if (standard.ized) {
    point.btw = ah.analysis[ah.analysis$x >= range(x$x)[1]*x.sd+x.mean & ah.analysis$x <= range(x$x)[2]*x.sd+x.mean,]
    point.out = ah.analysis[ah.analysis$x <= range(x$x)[1]*x.sd+x.mean | ah.analysis$x >  range(x$x)[2]*x.sd+x.mean,]
    
  } else {
    point.btw = ah.analysis[ah.analysis$x >= range(x$x)[1] & ah.analysis$x <= range(x$x)[2],]
    point.out = ah.analysis[ah.analysis$x <= range(x$x)[1] | ah.analysis$x > range(x$x)[2],]
  }
  
  
  # Plot the points that are from the GAM 
  if (!linear.only) {
    title=paste("Yr =",paste(my.eco.evo.df$yr1[i],my.eco.evo.df$yr2[i],collapse = "-"),
                "\n p-value of G =",round(sumglm$coefficients[3,4]/2,3), sep = ' ')
  } else {title = paste("Yr =",paste(my.eco.evo.df$yr1[i],my.eco.evo.df$yr2[i],collapse = "-"), sep = ' ')}
  plot(point.out$y~point.out$x,
       pch =21, 
       col = "black", 
       bg = "black",
       xlab = "PC1",
       ylab = "Survival",
       xlim = range(ah.analysis$x),
       ylim=c(0,1),
       main = title,
       cex = point.out$app.surv+.1)
  
  
  # Add text that was chosen as the boundary 
  if (standard.ized) {
    text(x = range(x$x*x.sd+x.mean)[1], y = .9,labels = round(range(x$x*x.sd+x.mean),2)[1], pos = 4)
    text(x = range(x$x*x.sd+x.mean)[2], y = .9,labels = round(range(x$x*x.sd+x.mean),2)[2], pos = 2)
    
  } else {
    text(x = range(x$x)[1],y = .9,labels = round(range(x$x),2)[1], pos = 4)
    text(x = range(x$x)[2],y = .9,labels = round(range(x$x),2)[2], pos = 2)
  }
  
  # Add line of GAM (fitness function)
  lines(ah.analysis$x, 
        ah.analysis$y)
  
  # Plot points that are chosen between the 2 fitness max 
  points(x = point.btw$x,
         y = point.btw$y, 
         cex = point.btw$app.surv+.1, 
         pch =21, col = "red", bg = "red")
  
  dpc1 = density(ah.analysis$x)
  # This is the range of x values (maximum of the atual data that I have)
  abline(v = range(x$x), 
         lty = 2, col = "blue")
  # This is the actual peak computed 
  # abline(v = c(my.eco.evo.df$local.max.peak1[i],
  #              my.eco.evo.df$local.max.peak2[i]), col = "blue", lty = 3, lwd = 3)
  
  newx2 <- seq(min(ah.analysis$x), 
               max(ah.analysis$x), 
               length.out = 1000)
  z2 <- predict(model.list[[i]], 
                newdata=list(x = newx2), 
                se.fit = TRUE)
  invlogit <- function(x){exp(x)/(exp(x) + 1)}
  yhat <- invlogit(z2$fit)
  upper <- invlogit(z2$fit + z2$se.fit)
  lower <- invlogit(z2$fit - z2$se.fit)
  mfakedata = data.frame(newx2,
                         z2,
                         yhat,
                         upper,
                         lower)
  # Add error to the GAM 
  lines(newx2, upper, lty = 2)
  lines(newx2, lower, lty = 2)
  xx <- c(newx2, rev(newx2))
  yy <- c(c(upper),
          rev(c(lower)))
  
  # Add polygon to make it cute 
  polygon(xx, yy,
          col = adjustcolor( "grey", alpha.f = 0.2),
          border = NA)
  
  # Add the line that correspond to the logistic regression 
  if (standard.ized) {
    newx = data.frame(x = seq(min(x$x*x.sd+x.mean),
                              max(x$x*x.sd+x.mean),
                              length.out = 100))
    newx$x2 = newx$x^2
    newx$y = predict(model.st,
                     newx, 
                     type="response") 
    logi = function(x){exp(x)/(1+exp(x))}
    newx$y.logi = logi(newx$y)
    # Adding the center part of the logistic regression 
    lines(x = newx$x,
          y = newx$y, 
          col = "red",
          lwd = 2, 
          ylim = c(0,1))
    newx2 = data.frame(x = seq(min(ah.analysis$x),
                               max(ah.analysis$x),
                               length.out = 100))
    
  } else {
    newx = data.frame(x = seq(min(x$x),
                              max(x$x),
                              length.out = 100))
    newx$x2 = newx$x^2
    newx$y = predict(model,
                     newx, 
                     type="response") 
    logi = function(x){exp(x)/(1+exp(x))}
    newx$y.logi = logi(newx$y)
    
    # Adding the center part of the logistic regression 
    lines(x = newx$x,
          y = newx$y, 
          col = "red",
          lwd = 2, 
          ylim = c(0,1))
    newx2 = data.frame(x = seq(min(ah.analysis$x),
                               max(ah.analysis$x),
                               length.out = 100))
  }
  
  min.fit.trait = newx$x[which(newx$y == min(newx$y))]
  min.fit = min(newx$y)
  list.min.fit = c(list.min.fit,min.fit)
  list.min.fit.trait = c(list.min.fit.trait,min.fit.trait)
  
  if (standard.ized) {
    newx2.1 = newx2[newx2$x <= min(x$x*x.sd+x.mean),,drop=FALSE]
    newx2.2 = newx2[newx2$x >= max(x$x*x.sd+x.mean),,drop=FALSE]
  } else {
    newx2.1 = newx2[newx2$x <= min(x$x),,drop=FALSE]
    newx2.2 = newx2[newx2$x >= max(x$x),,drop=FALSE]
  }
  
  # Adding both part (away from the center) of the logistic regression 
  if (standard.ized) {
    newx2.1$x2 = newx2.1$x^2
    newx2.2$x2 = newx2.2$x^2
    newx2.1$y = predict(model.st,
                        newx2.1, 
                        type="response") 
    newx2.2$y = predict(model.st,
                        newx2.2, 
                        type="response") 
    newx2.1$y.logi = logi(newx2.1$y)
    newx2.2$y.logi = logi(newx2.2$y)
    
    lines(x = newx2.1$x,
          y = newx2.1$y, 
          col = "red",
          lty = 3,
          lwd=2, ylim = c(0,1))
    lines(x = newx2.2$x,
          y = newx2.2$y, 
          col = "red",
          lty = 3,
          lwd=2, ylim = c(0,1))
  } else {
    newx2.1$x2 = newx2.1$x^2
    newx2.2$x2 = newx2.2$x^2
    newx2.1$y = predict(model,
                        newx2.1, 
                        type="response") 
    newx2.2$y = predict(model,
                        newx2.2, 
                        type="response") 
    newx2.1$y.logi = logi(newx2.1$y)
    newx2.2$y.logi = logi(newx2.2$y)
    
    lines(x = newx2.1$x,
          y = newx2.1$y, 
          col = "red",
          lty = 3,
          lwd=2, ylim = c(0,1))
    lines(x = newx2.2$x,
          y = newx2.2$y, 
          col = "red",
          lty = 3,
          lwd=2, ylim = c(0,1))
    
  }
  
  
  
  
  
  # Add density distribution 
  par(new=TRUE)
  plot(dpc1, #axes=FALSE, ann=FALSE,
       zero.line=FALSE, col = "blue", 
       cex = .6, 
       main = "", 
       xlim = range(ah.analysis$x),
       axes=FALSE, 
       ann=FALSE, frame =FALSE)
  axis(side=4, at = pretty(range(dpc1$y)))
  
  
  # my.line=function(x, model) {
  #   out= model$coefficients[1] +
  #     c(model$coefficients[2]*x) +
  #     c(model$coefficients[3]*x^2)
  #   logi=function(x){exp(x)/(1+exp(x))}
  #   return(
  #     logi(out))
  # }
  
  
  # curve(logi, newx, 
  #       col = "red",
  #       lwd=2,
  #       add=TRUE,
  #       xlim = range(ah.analysis$x))
  # lines(x = x,y = my.line(x,model), col = "red",lwd=2)
  # curve(predict(model,data.frame(x=x)), col = "red",lwd=2,add=TRUE)
  
  library(aod)
  if (linear.only) {
    wald.test(b = coef(object=model), 
              Sigma = vcov(object=model), 
              Terms = 2)
  } else{
    wald.test(b = coef(object=model), 
              Sigma = vcov(object=model), 
              Terms = 3)
    
  }
  anova(object=model, test="Chisq")
  summary(model)
  
  coef(summary(model))
  confint.default(model)
  confint(model)
} # End of for(i in 1:nrow(my.eco.evo.df)){ 

library(data.table)
dt<-rbindlist(fit.grad.table)
dt = as.data.frame(dt)
if (linear.only) {
  row.names(dt) <- c("x_04_05", 
                     "x_05_06", 
                     "x_06_07",
                     "x_07_08",   
                     "x_08_09",  
                     "x_09_10", 
                     "x_10_11")
  
} else {
  row.names(dt) <- c("x_04_05", "x2_04_05", 
                     "x_05_06", "x2_05_06", 
                     "x_06_07", "x2_06_07",
                     "x_07_08", "x2_07_08",   
                     "x_08_09", "x2_08_09",  
                     "x_09_10", "x2_09_10", 
                     "x_10_11", "x2_10_11")
  
}
# dput(paste("x",substr(my.eco.evo.df$yr1,3,4),substr(my.eco.evo.df$yr2,3,4), sep = "_"))
# dput(paste("x2",substr(my.eco.evo.df$yr1,3,4),substr(my.eco.evo.df$yr2,3,4), sep = "_"))
dt$lin.p.1.t = dt$lin.p/2 
dt$lgt.p.1.t = dt$lgt.p/2 

# Look at only significant gradients for a one tail p value (greater than 0, because a>0 is shaped like a bowl)
dt[which(dt$lgt.p.1.t<=0.05),]
dt.sign = dt[which(dt$lgt.p.1.t<=0.05),]
row.names(model.list.logistic) <- c("x_04_05", 
                                    "x_05_06", 
                                    "x_06_07",
                                    "x_07_08",   
                                    "x_08_09",  
                                    "x_09_10", 
                                    "x_10_11")
model.list.logistic = as.data.frame(model.list.logistic)
model.list.logistic$sample.n.0 = do.call(rbind, sample.size)[,1] # alternative 1. Using R base
model.list.logistic$sample.n.1 = do.call(rbind, sample.size)[,2] # alternative 1. Using R base
model.list.logistic$sample.n.tot = model.list.logistic$sample.n.0 + model.list.logistic$sample.n.1

df.final.df.new.analysis = do.call(rbind,final.df.new.analysis)

# Save everything else ------------------------------------------------------

save(model.list.logistic, # Summary of all the GLM models (logistic) 
     glm.model.list, # This is the model output of all the GLMs calculated from the logistic regression 
     dt, # Ran from Janzen program to calculate selection gradients  
     dt.sign, # Only the significant gradients 
     pseudr, # pseudo r squared of the logistic models 
     list.min.fit, 
     list.min.fit.trait, 
     gofitfit, # Goodness of fit of the logistic models using Hosmer–Lemeshow test
     standard.ized = standard.ized,
     linear.only = linear.only,
     orthogonal.x = orthogonal.x,
     original.x = original.x,
     final.df.new.analysis,
     df.final.df.new.analysis, 
     file = "~/Dropbox/finch_recap/v2/data/Selection coefficients/logistic_gam_models.selection.RData") # Where is it! 

load("~/Dropbox/finch_recap/v2/data/Selection coefficients/logistic_gam_models.selection.RData", verbose = TRUE)
# Model GLM with all data -------------------------------------------------
load("~/Dropbox/finch_recap/v2/data/creating_data/climate.year.month.database.RData", 
     verbose = TRUE) # created from ~/Dropbox/finch_recap/v2/data/creating_data/Climate_all_time_Santa_cruz.R
db.climate.year # Per year precipitation, sum, mean, var, sd, se 

as.factor(df.final.df.new.analysis$y)

# This is to find the rows in each year to get the rain data  
row.precip = match(c(as.numeric(as.character(df.final.df.new.analysis$yea.var)) -2),db.climate.year$Year.marco)
# I decided to look for the rain data in the complete year (wet and dry combined) 
df.final.df.new.analysis$precipitation.yr.before = db.climate.year[row.precip,"sum.prec"]
# I decided to look for the rain data in the dry season only 
df.final.df.new.analysis$prcp.yr.before.dry = dry.rain.our.study[row.precip,"sum.prec"]
df.final.df.new.analysis$prcp.yr.before.wet = wet.rain.our.study[row.precip,"sum.prec"]

df.final.df.new.analysis = df.final.df.new.analysis[df.final.df.new.analysis$yea.var != 2007,]
plot(y~my.x, 
     data = df.final.df.new.analysis)
df.final.df.new.analysis$x = df.final.df.new.analysis$my.x
df.final.df.new.analysis$x2 = df.final.df.new.analysis$my.x^2
df.final.df.new.analysis$yea.var  = as.factor(df.final.df.new.analysis$yea.var)
glm.out1 =  glm(y~ x,                  data = df.final.df.new.analysis,family = "binomial")
glm.out2 =  glm(y~ x + x2,             data = df.final.df.new.analysis,family = "binomial")
glm.out3 =  glm(y~ x + x2 + yea.var,   data = df.final.df.new.analysis,family = "binomial")
glm.out4 =  glm(y~ x + x2 * yea.var,   data = df.final.df.new.analysis,family = "binomial")
glm.out5 =  glm(y~ x + x2 : yea.var,   data = df.final.df.new.analysis,family = "binomial")
glm.out6 =  glm(y~(x + x2) : yea.var,  data = df.final.df.new.analysis,family = "binomial")
glm.out7 =  glm(y~ x * yea.var + x2 * yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out8 =  glm(y~ x:x2:yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out9 = glm(y~ x:yea.var+ x2:yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out10 = glm(y~ x:yea.var+ x2:yea.var + yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out11 =  glm(y~(x + x2) * yea.var,  data = df.final.df.new.analysis,family = "binomial")
glm.out12 =  glm(y~ x*x2*yea.var, data = df.final.df.new.analysis,family = "binomial")
glm.out13 =  glm(y~ x*x2*yea.var*precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out14 =  glm(y~ x + x2 * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out15 =  glm(y~ x + x2 : precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out16 =  glm(y~ (x + x2)* precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
# glm.out3 =  glm(y~ x + x2,             data = df.final.df.new.analysis,family = "binomial")
summary(glm.out11)
summary(glm.out16)

glm.out17 =  glm(y~ x + x2 + precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out18 =  glm(y~ x + x2 * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out19 =  glm(y~ x + x2 : precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out20 =  glm(y~(x + x2) : precipitation.yr.before,  data = df.final.df.new.analysis,family = "binomial")
glm.out21 =  glm(y~ x * precipitation.yr.before + x2 * precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out22 =  glm(y~ x:x2:precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out23 = glm(y~ x:precipitation.yr.before+ x2:precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out24 = glm(y~ x:precipitation.yr.before+ x2:precipitation.yr.before + precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out25 =  glm(y~(x + x2) * precipitation.yr.before,  data = df.final.df.new.analysis,family = "binomial")
glm.out26 =  glm(y~ x*x2*precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out27 =  glm(y~ x*x2*precipitation.yr.before*precipitation.yr.before, data = df.final.df.new.analysis,family = "binomial")
glm.out28 =  glm(y~ x + x2 * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out29 =  glm(y~ x + x2 : precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
table(df.final.df.new.analysis$yea.var)
glm.out30 =  glm(y~ (x + x2)* precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")

glm.out31 =  glm(y~ (x + x2)* prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out32 =  glm(y~ x + (x2)* prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out33 =  glm(y~ x + x2 + prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out34 =  glm(y~ x + x2: prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")


glm.out28 =  glm(y ~  x + x2  * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out32 =  glm(y ~  x + (x2)* prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out35.1 =  glm(y~ x + (x2)* prcp.yr.before.wet,   data = df.final.df.new.analysis,family = "binomial")
summary(glm.out28)
summary(glm.out32)
summary(glm.out35.1)


summary(glm.out31)
summary(glm.out32)
summary(glm.out33)
summary(glm.out34)
summary(glm.out35.1)


df.final.df.new.analysis$prc.bfr = df.final.df.new.analysis$precipitation.yr.before
glm.out102=summary(glm(y~x2 * prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out103=summary(glm(y~x + x2 * prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out104=summary(glm(y~(x + x2) * prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out105=summary(glm(y~ x2+ prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out106=summary(glm(y~ x + x2,data=df.final.df.new.analysis,family="binomial"))
glm.out107=summary(glm(y~ x2,data=df.final.df.new.analysis,family="binomial"))
glm.out108=summary(glm(y~ x+x2+prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out109=summary(glm(y~ x*prc.bfr+x2,data=df.final.df.new.analysis,family="binomial"))
glm.out110=summary(glm(y~ prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out111=summary(glm(y~ 1,data=df.final.df.new.analysis,family="binomial"))
glm.out112=summary(glm(y~ x*prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out113=summary(glm(y~ x+prc.bfr,data=df.final.df.new.analysis,family="binomial"))
glm.out114=summary(glm(y~ x,data=df.final.df.new.analysis,family="binomial"))

mod.l=list(glm.out102,glm.out103,glm.out104,
           glm.out105,glm.out106,glm.out107,glm.out108,
           glm.out109,glm.out110,glm.out111,
           glm.out112,glm.out113,glm.out114)
numb.modl = unlist(lapply(lapply(mod.l,coefficients),nrow))

test.all.var = do.call("rbind", (lapply(mod.l,coefficients)))
rown= row.names(test.all.var)

test.all.var = data.frame(test.all.var)
test.all.var$var= rown
mlistt = NULL
for (i in 1:length(numb.modl)) {
  mlistt = c(mlistt,rep(i,numb.modl[i]))
}
test.all.var$model = (mlistt)
write.csv(test.all.var,"~/Desktop/New article to merge LUKE/Images/Tables/glm.model.p.values.coefficients.csv",na = "")

# model with just interaction 
# model GLMM precipitation (fixed), year (random)
summary(glm.out1)
summary(glm.out2)
summary(glm.out3)
summary(glm.out4) 
summary(glm.out5) 
summary(glm.out6) 
summary(glm.out7) 
summary(glm.out8) 
summary(glm.out9) 
summary(glm.out10) 
summary(glm.out11) 
summary(glm.out12) 
summary(glm.out13) 
summary(glm.out14) 
summary(glm.out15) 
summary(glm.out16) 

summary(glm.out17)
summary(glm.out18)
# summary(glm.out19)
# summary(glm.out20) 
summary(glm.out21) 
# summary(glm.out22) 
# summary(glm.out23)
# summary(glm.out24) 
summary(glm.out25) 
# summary(glm.out26) 
# summary(glm.out27) 
summary(glm.out28) 
# summary(glm.out29) 
summary(glm.out30) 


library(MuMIn)
old.options=options()
options(na.action = na.fail)
# see https://rcompanion.org/rcompanion/d_04.html 
options(contrasts = c("contr.sum", "contr.poly"))
### needed for type III tests

glm.out28 =  glm(y ~  (x + x2) * precipitation.yr.before,   data = df.final.df.new.analysis,family = "binomial")
glm.out32 =  glm(y ~  (x + x2) * prcp.yr.before.dry,   data = df.final.df.new.analysis,family = "binomial")
glm.out35.1 =  glm(y~ (x + x2) * prcp.yr.before.wet,   data = df.final.df.new.analysis,family = "binomial")
model.sel28 = dredge(glm.out28)
model.sel32 = dredge(glm.out32)
model.sel35 = dredge(glm.out35.1)
# Automated model selection
model.sel0 = dredge(glm.out7)
model.sel1 = dredge(glm.out30)
model.sel2 = dredge(glm.out31)
model.sel2 = dredge(glm.out31)
# model.sel2 = dredge(glm.out13) #too much to look at! 
model.avg(model.sel0,subset = delta < 4)
model.avg(model.sel,subset = delta < 4)
att.mod0 =attr(model.sel0,"model.calls")
att.mod1 =attr(model.sel1,"model.calls")
att.mod2 =attr(model.sel2,"model.calls")
mod.0 = gsub(",     ", replacement = "", paste(format(att.mod0)))
mod.0 = gsub('\"', replacement = "", mod.0)
mod.1 = gsub(",     ", replacement = "", paste(format(att.mod1)))
mod.1 = gsub('\"', replacement = "", mod.1)
mod.2 = gsub(",     ", replacement = "", paste(format(att.mod2)))
mod.2 = gsub('\"', replacement = "", mod.2)

model.sel0$models = mod.0
model.sel1$models = mod.1
model.sel2$models = mod.2

write.csv(model.sel0,"~/Desktop/New article to merge LUKE/Images/Tables/glm.all.data_year.csv",na = "")
write.csv(model.sel1,"~/Desktop/New article to merge LUKE/Images/Tables/glm.all.data_precip.csv",na = "")

write.csv(model.sel2,"~/Desktop/New article to merge LUKE/Images/Tables/glm.all.data_year_precip.csv")
model.sel1
model.sel2
str(model.sel1)
coeffs(model.sel1)
attr(model.sel1,"coefTable")
row.names(model.sel0[model.sel0$delta <=4,])
attr(model.sel0,"coefTable")
attr(model.sel0,"model.calls")$`8`
attr(model.sel0,"model.calls")$`4`
attr(model.sel0,"model.calls")$`3`
attr(model.sel0,"model.calls")$`16`
attr(model.sel0,"model.calls")$`24`

attr(model.sel1,"model.calls")$`22`
summary(glm(formula = y ~ x + x2 + yea.var + 1, family = "binomial", 
            data = df.final.df.new.analysis))
summary(glm(formula = y ~ x + x2 + 1, family = "binomial", data = df.final.df.new.analysis))
summary(glm(formula = y ~ x2 + 1, family = "binomial", data = df.final.df.new.analysis))
summary(glm(formula = y ~ x + x2 + yea.var + x:yea.var + 1, family = "binomial", 
            data = df.final.df.new.analysis))
summary(glm(formula = y ~ x + x2 + yea.var + x2:yea.var + 1, family = "binomial", 
            data = df.final.df.new.analysis))

library(car)
Anova(glm.out1,type = "III", test.statistic = "Wald")
Anova(glm.out2,type = "III", test.statistic = "Wald")
Anova(glm.out3,type = "III", test.statistic = "Wald")
Anova(glm.out4,type = "III", test.statistic = "Wald")
Anova(glm.out5,type = "III", test.statistic = "Wald")
Anova(glm.out6,type = "III", test.statistic = "Wald")
Anova(glm.out7,type = "III", test.statistic = "Wald")
Anova(glm.out8,type = "III", test.statistic = "Wald")
Anova(glm.out9,type = "III", test.statistic = "Wald")
Anova(glm.out10,type = "II", test.statistic = "Wald")
Anova(glm.out10,type = "III", test.statistic = "Wald")
Anova(glm.out11,type = "III", test.statistic = "Wald")
Anova(glm.out12,type = "III", test.statistic = "Wald")
Anova(glm.out13,type = "III", test.statistic = "Wald")
Anova(glm.out14,type = "III", test.statistic = "Wald")
Anova(glm.out15,type = "III", test.statistic = "Wald")
Anova(glm.out16,type = "III", test.statistic = "Wald")

Anova(glm.out17,type = "III", test.statistic = "Wald")
Anova(glm.out18,type = "III", test.statistic = "Wald")
Anova(glm.out19,type = "III", test.statistic = "Wald")
Anova(glm.out20,type = "III", test.statistic = "Wald")
Anova(glm.out21,type = "III", test.statistic = "Wald")
Anova(glm.out22,type = "III", test.statistic = "Wald")
Anova(glm.out23,type = "III", test.statistic = "Wald")
Anova(glm.out24,type = "III", test.statistic = "Wald")
Anova(glm.out25,type = "III", test.statistic = "Wald") #***
Anova(glm.out26,type = "II", test.statistic = "Wald")
Anova(glm.out26,type = "III", test.statistic = "Wald")
Anova(glm.out27,type = "III", test.statistic = "Wald")
Anova(glm.out28,type = "III", test.statistic = "Wald")
Anova(glm.out29,type = "III", test.statistic = "Wald")
Anova(glm.out30,type = "III", test.statistic = "Wald")

attr(model.sel1,"model.calls")$`32`

modmod =glm(formula = y ~ precipitation.yr.before + x2 + precipitation.yr.before:x2 + 
              1, family = "binomial", data = df.final.df.new.analysis)
Anova(modmod, type = "III", test.statistic = "Wald")

modmod =glm(formula = y ~ precipitation.yr.before + x + x2 + precipitation.yr.before:x2 + 
              1, family = "binomial", data = df.final.df.new.analysis)
Anova(modmod, type = "III", test.statistic = "Wald")

modmod =glm(formula = y ~ precipitation.yr.before + x + x2 + precipitation.yr.before:x + 
              precipitation.yr.before:x2 + 1, family = "binomial", data = df.final.df.new.analysis)
Anova(modmod, type = "III", test.statistic = "Wald")

Anova(glm.out1,type = "II", test.statistic = "Wald")
Anova(glm.out2,type = "II", test.statistic = "Wald")
Anova(glm.out3,type = "II", test.statistic = "Wald")
Anova(glm.out4,type = "II", test.statistic = "Wald")
Anova(glm.out5,type = "II", test.statistic = "Wald")
Anova(glm.out6,type = "II", test.statistic = "Wald")
Anova(glm.out7,type = "II", test.statistic = "Wald")
Anova(glm.out8,type = "II", test.statistic = "Wald")
Anova(glm.out9,type = "II", test.statistic = "Wald")
Anova(glm.out10,type = "II", test.statistic = "Wald")
Anova(glm.out11,type = "II", test.statistic = "Wald")
Anova(glm.out12,type = "II", test.statistic = "Wald")

anova(glm.out1,glm.out2,
      glm.out3,glm.out4,glm.out5,glm.out5,glm.out6,glm.out7,glm.out8,glm.out9,glm.out10)


library(lme4)
df.final.df.new.analysis$pcr.scale = scale(df.final.df.new.analysis$prc.bfr)
rd.modl = glmer(formula = y ~ prc.bfr + 
                  x + x2 + 
                  prc.bfr:x + prc.bfr:x2 + 
                  (1 | my.band),
                family = "binomial", data = df.final.df.new.analysis)
rd.modl1 = glmer(formula = y ~ x + x2 + prc.bfr+
                   (1 | my.band), family = "binomial", 
                 data = df.final.df.new.analysis)
rd.modl1 = glmer(formula = y ~ x + x2* prc.bfr+
                   (1 | my.band), family = "binomial", 
                 data = df.final.df.new.analysis)
summary(rd.modl)
summary(rd.modl1)
dredge(rd.modl) # Takes a long time 
library(sjPlot)
sjp.glmer(rd.modl, type = "fe")

par(mfrow=c(1,1))

# plot(glm.out)
newx <- seq(from = min(df.final.df.new.analysis$my.x), 
            to = max(df.final.df.new.analysis$my.x), 
            length.out = 2000)

z1 <- predict(glm.out1, 
              newdata = data.frame(x = newx,  
                                   x2 = newx^2, 
                                   yea.var = 2009), 
              se.fit = TRUE)
# Function that will transforme the Y values from linear scale to [0,1]
invlogit <- function(x){exp(x)/(exp(x) + 1)}
# This is the actual transformed data 
yhat <-  invlogit(z1$fit)
upper <- invlogit(z1$fit + z1$se.fit)
lower <- invlogit(z1$fit - z1$se.fit)

# satisfied...= "n"
# while(satisfied... != "y"){
# Here is the plot of the fitness function using the evenly spaced data (newx) and the response to it (yhat)
plot(newx, yhat, type="l", 
     ylim = c(0,1), 
     xlab = "PC1",
     main = "GLM all data")
# Adding error 
lines(newx, upper, lty = 2)
lines(newx, lower, lty = 2)
points(df.final.df.new.analysis$my.x,df.final.df.new.analysis$y)
# par(new = TRUE)



# LDA ---------------------------------------------------------------------
library(MASS)
library(klaR)
mld = lda(y~my.mbl +my.mbw + my.mbd,data = df.final.df.new.analysis)
mqd = qda(y~my.mbl +my.mbw + my.mbd,data = df.final.df.new.analysis)
mld
mqd
plot(mld)
# par(mfrow=c(4,2))
# for (j  in 1:length(unique(df.final.df.new.analysis$yea.var))) {
#   all.y = unique(df.final.df.new.analysis$yea.var)
#   my.df.for = df.final.df.new.analysis[df.final.df.new.analysis$yea.var == all.y[j],]
# partimat(as.factor(my.df.for$y) ~ my.mbl +my.mbw + my.mbd,
#          data = my.df.for, method = "lda", plot.matrix = TRUE)
# }
if (pdf.output) {
  pdf("~/Dropbox/finch_recap/v2/figures/quadratic_discriminant_function.pdf")
}
# partimat(as.factor(df.final.df.new.analysis$y) ~ my.mbl +my.mbw + my.mbd,
#          data = df.final.df.new.analysis, method = "qda", prec = 50, 
#          col.correct = "black",
#          col.wrong = "red", 
#          imageplot = TRUE,
#          # gs = as.numeric((df.final.df.new.analysis$y)),
#          pch.mean = 1:length(levels(as.factor(df.final.df.new.analysis$y))),
#          cex.mean = 1.4,
#          # image.colors = couleur,
#          col.mean = "yellow",
#          plot.matrix = TRUE)
if (pdf.output) {
  dev.off()
}


# Relation with environment -----------------------------------------------
my.eco.evo.df$preci.yr1

# Extract only linear 
if (!linear.only) {
  lin = dt[seq(1,nrow(dt),by = 2),]
} else {lin = dt}

# lin = lin[lin$lgt.p.1.t<= 0.05,]
plot(lin$bavggrd~my.eco.evo.df$preci.yr2)
plot(lin$bavggrd.mean.st~my.eco.evo.df$preci.yr2)

# Extract only quadratic 
quad = dt[seq(2,nrow(dt),by = 2),]

# Here I wanted to know what was the sample size for males and females all at once. 
yr = c(2004:2006,2007, 2008:2011)
source('src/initialize.R')
for(i in 1:c(length(yr)-1)){
  yr.list <- c(yr[i],yr[i+jump])
  mdat <- prep.data(sp.keep = sp.list,
                    center = FALSE,yr.keep = yr.list,
                    site.keep = site.list,flip.pc1 = TRUE,
                    adults.only = FALSE,gr.sel = NULL,
                    valley.fortis = FALSE, valley.fortis.fuli = FALSE,
                    peaks.fortis.fuli = FALSE,two.univariate.traits = FALSE,
                    keep.resid = FALSE, splines = splines,
                    keep.last.year.data = FALSE,gam.analysis = TRUE,
                    K.nots = 5, # Manually select the number of knots  
                    traits = c('PC1'))#, 'PC2')) # This is only for the spline 
  tab.int = data.frame(surv = mdat$X[,2],
                       sex = droplevels(mdat$ind.vars$age),
                       pc1 = mdat$ind.vars$pc1)
  yearly.number.of.id = c(yearly.number.of.id,list(tab.int))
} 
new.yearly.list = NULL
for (i in 1:7) {
  vector.find.cut.each.year = yearly.number.of.id[[i]]$pc1>= all.ranges.x[[i]][1]  & 
    yearly.number.of.id[[i]]$pc1<= all.ranges.x[[i]][2]
  new.yearly.list = c(new.yearly.list,list(yearly.number.of.id[[i]] [vector.find.cut.each.year,]))
}

# Without cutting anything from the data 
lapply(yearly.number.of.id, function(x) table(x$sex))
lapply(yearly.number.of.id, function(x) table(x$surv))
lapply(new.yearly.list, function(x) table(x$sex))
lapply(new.yearly.list, function(x) table(x$surv))

lapply(yearly.number.of.id, nrow)
par(mfrow=c(2,2))
if (!linear.only) {
  lm.o1 = lm(quad$bavggrd~my.eco.evo.df$preci.yr1)
  summary(lm.o1)
  plot(quad$bavggrd~my.eco.evo.df$preci.yr1);abline(lm.o1)
  
  lm.o2 = lm(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr1)
  summary(lm.o2)
  plot(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr1);abline(lm.o2)
  
  lm.o3 = lm(quad$bavggrd~my.eco.evo.df$preci.yr2)
  summary(lm.o3)
  plot(quad$bavggrd~my.eco.evo.df$preci.yr2);abline(lm.o3)
  
  lm.o4 = lm(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr2)
  summary(lm.o4)
  plot(quad$bavggrd.mean.st~my.eco.evo.df$preci.yr2);abline(lm.o4)
}
if(pdf.output){
  dev.off()
  dev.off()
  dev.off()
  dev.off()
}







