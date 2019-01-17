# Load the data  ----------------------------------------------------------
load("./output/biotic.factors.on.survival_Andrew_meeting_changing_PCA_SCORE_for_only_fortis.RData", 
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

bartlett.test(pc1 ~ gr,data=compar)
# If not normally distributed 
car::leveneTest(pc1 ~ gr,data=compar)

# plot(mixturemodel.out.groups$pc1, col = mixturemodel.out.groups$gr)
par(mfrow=c(3,4))
for (i in 1:12) {
  plot(ores.beak[[i]],whichplots = 2, main2 = c(2005:2018)[i])
}
par(mfrow=c(1,1))

par(mfrow= c(2,1))
mydf = data.frame(y = newlist[[1]], x = newx)
where.peak1 = which(mydf$y == max(mydf$y[mydf$x>my.eco.evo.df$local.max.peak1[1]-.01  & 
                                           mydf$x<my.eco.evo.df$local.max.peak1[1]+.1]))

sub.newx = newx[where.peak1:nrow(mydf)]
x = sub.newx
q.out = lm(mydf$y[where.peak1:nrow(mydf)]~poly(x,2))
summary(q.out)

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
  original.x = c(original.x,list(sub$x))
  my.x = sub$x
  my.mbl = sub$mbl
  my.mbw = sub$mbw
  my.mbd = sub$mbd
  my.band = sub$band
  y = sub$app.surv
  yea.var=rep(my.eco.evo.df[i,"yr2"],length(y))
  new.analysis.all.years = data.frame(my.x,y,yea.var,my.mbl,my.mbw,my.mbd,my.band)
  final.df.new.analysis = c(final.df.new.analysis,list(new.analysis.all.years))
  
  # Preparing the data X (quadratic or not, orthongonal or not)
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
  
  # According to endler (p.255), if the covariation between X and X^2 is 0, 
  # there is no skweness and therefore the univariate regression and 
  # multivariate regression coefficients are the same 
  kurtosis(x$x)
  skewness(x$x) 
  cov(x$x,(x$x-mean(x$x))^2)
  
  galtonskew.proc(x$x)
  
  sample.size = c(sample.size,
                  list(table(y)))
  ch = cbind(x= rep(1,length(y)),y)
  freq = rep(1,length(y))
  
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
  gofit = HLtest(model)
  cat("The p-value of the GOF is:",round(gofit$p.value,2),"\n")
  
  pseudr = c(pseudr,pseudoR2)
  gofitfit = c(gofitfit,list(gofit))
  
  
  model.list.logistic = rbind(model.list.logistic,
                              as.table(coefficients(model)))
  sumglm = summary(model)
  
  
  # Compute Gradients -----------------------------------------
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
     file = "./output/logistic_gam_models.selection.RData") # Where is it! 

