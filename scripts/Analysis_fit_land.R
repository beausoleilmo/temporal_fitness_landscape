# # # # # # # # # # # # # # # # # # # # 
# Author: Marc-Olivier Beausoleil 
# Date: January 17, 2019
# McGill University 
# Analysis of the dynamics of Darwin's finches' fitness landscapes  
# # # # # # # # # # # # # # # # # # # # 

# Preparation of variables and data  --------------------------------------
source('./scripts//initialize.R')
source("./scripts/logit.r")
source('./scripts/PCA_custom.R')
load('./data/bird.data.RData', verbose=TRUE)

# Select the variables to run the scripts ---------------------------------
save.data = "./output/biotic.factors.on.survival_Andrew_meeting_changing_PCA_SCORE_for_only_fortis.RData"

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
exp.lambda = exp(mean(c(-4,-5,-4,
                        -4,-4,-4,
                        -4,-4,-3,
                        -6,-13,0)))

if(pdf.output){
  pdf("./output/my.pca.just.fortis.pdf",height = 6,width = 15)
  par(mfrow=c(1,3))
  res.pca.for2 = vegan::rda(bird.data[bird.data$Species1=="fortis", 
                                      c("MedianBeakLength", "MedianBeakWidth", 
                                        "MedianBeakDepth")])
  res.pca.for3 = vegan::rda(bird.data[bird.data$Species1=="fortis" & bird.data$Site=="El Garrapatero",
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
  
  tab.sp.eg= table(bird.data[bird.data$Site=="El Garrapatero","Species1"], 
                   bird.data[bird.data$Site=="El Garrapatero","Year"])
  write.csv(tab.sp.eg,"./output/nb.sp.eg.csv")
  
  ores = mixtools::normalmixEM(bird.data[bird.data$Species1=="fortis","PC1"],
                               sigma = NULL, 
                               mean.constr = NULL, sd.constr = NULL,
                               epsilon = 1e-15, maxit = 1000, maxrestarts=50, 
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
  
  # res.pcabs is found in load('./data/bird.data.RData', verbose=TRUE)
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
  pdf(paste0("./output/my.fit.land.short2_jump",
             jump,"_",gsub('([[:punct:]])|\\s+','_',site.list),".pdf"),
      height = 7,
      width = 8)
}

# Function finding maximum and minimum of the fitness function  -----------
par(mfrow=c(1,1))
# For loop that will calculate the GAM, show the landscape and 
# will let you select what is the maximum and minimum of the function 
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
                      gr.sel = NULL,# "small" or "big", this is to create 2 groups in the fortis population. # applicable for fortis only. This will let you select "big" or "small" morph
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
                      K.nots = 5, # Manually select the number of knots (for spline)
                      traits = c('PC1')) # This is only for the spline 

    # This is to calcualte the GAM 
    # getting response variable and the explanatory variable 
    y = as.vector(mdat$X[,2])
    x = c(mdat$ind.vars$pc1) # Works for everything execpt 2007-2008
    mbd = c(mdat$ind.vars$mbd) 
    mbl = c(mdat$ind.vars$mbl) 
    mbw = c(mdat$ind.vars$mbw) 
    band = as.character(mdat$ind.vars$band) 
    year.var = rep(yr.list[2],length(mdat$ind.vars$pc1))
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
    
    yr1 = substr(i + 2003, 3, 4)
    yr2 = substr(i + 2003 + jump, 3, 4)
    oldxlist = c(oldxlist, list(x))
    oldzlist = c(oldzlist, list(z$fitted.values))
    old_beak_L_list = c(old_beak_L_list, list(mbd))
    old_beak_W_list = c(old_beak_W_list, list(mbl))
    old_beak_D_list = c(old_beak_D_list, list(mbw))
    old_band_list = c(old_band_list, list(band))
    
    survived.list = c(survived.list, list(y))
    
    lambb.z = round(log(exp.lambda))
    
    # Getting new x that is spaced evenly respecting the GAM function (fitness function) 
    newx <- seq(from = min(mydata$x), 
                to = max(mydata$x), 
                length.out = 2000)
    
    # Using the model to generate the new response variable 
    z1 <- predict(z, 
                  newdata=list(x = newx), 
                  se.fit = TRUE)
    
    # This is the actual transformed data 
    yhat <- invlogit(z1$fit)
    upper <- invlogit(z1$fit + z1$se.fit)
    lower <- invlogit(z1$fit - z1$se.fit)
    
    # Here is the plot of the fitness function using the evenly spaced 
    # data (newx) and the response to it (yhat)
    plot(newx, yhat, type="l", 
         ylim = c(0,1), 
         xlab = "PC1",
         main = bquote(atop("GAM±1SE",
                            lambda * " = "*.(lambb.z) * 
                              ",  yr = "*.(yr1) *"-"* .(yr2))))
    # Adding error 
    lines(newx, upper, lty = 2)
    lines(newx, lower, lty = 2)
    
    # find the *minimum* of the fitness function from the gam by clicking on BOTH sides of the highest visible peak and  the minimum value between the 2 peaks (valley) in the GAM  
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
    newlist = c(newlist, list(yhat))
    
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
  
  eff = structure(c(36, 140, 212, 120, 52, 56, 
                    132, 300, 128, 120, 128, 
                    104, 30.2565802161513, 81.4673085182818, 
                    81.4673085182818, 117.180771285717
  ), .Names = c("2003", "2004", "2005", 
                "2006", "2007", "2008", 
                "2009", "2010", "2011", 
                "2012", "2013", "2014", 
                "2015", "2016", "2017", "2018"))
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
} # End of if(find.peaks.and.valleys){


# Save finding peaks ------------------------------------------------------
if(find.peaks.and.valleys){
  save(my.eco.evo.df,
       lm.out1, lm.out2,
       newlist,
       newx,
       model.list,
       oldxlist, oldzlist,
       old_beak_L_list, old_beak_W_list, old_beak_D_list,
       old_band_list,
       survived.list,
       full.data,
       file = save.data)
}